library(atakrig)
library(raster)

## part 0: data preparation, discretize data ----
rpath <- system.file("extdata", package="atakrig")
aod3k <- raster(file.path(rpath, "MOD04_3K_A2017042.tif"))
aod10 <- raster(file.path(rpath, "MOD04_L2_A2017042.tif"))

aod3k.d <- discretizeRaster(aod3k, 1500)
aod10.d <- discretizeRaster(aod10, 1500)
grid.pred <- discretizeRaster(aod3k, 1500, type = "all")

aod3k.d$areaValues$value <- log(aod3k.d$areaValues$value)
aod10.d$areaValues$value <- log(aod10.d$areaValues$value)


## part 1: gstat ordinary Kriging ----
aod3k.pt <- aod3k.d$areaValues
vgm.ok <- autofitVgm(aod3k.pt, ngroup=12, rd=0.75, fig = TRUE)$model
unknown.pt <- grid.pred$areaValues
coordinates(aod3k.pt) = ~centx+centy
coordinates(unknown.pt) = ~centx+centy
pred.ok <- gstat::krige(value ~ 1, aod3k.pt, newdata = unknown.pt,
                        model = vgm.ok, nmax = 10)
pred.ok$var1.pred <- exp(pred.ok$var1.pred)
gridded(pred.ok) <- TRUE
pred.ok.r <- raster(pred.ok[1])


## part 2: point-scale variogram deconvolution ----
# point-scale variogram from combined AOD-3k and AOD-10
aod.combine <- rbindDiscreteArea(x = aod3k.d, y = aod10.d)
vgm.ok_combine <- deconvPointVgm(aod.combine, model="Exp", ngroup=12, rd=0.75)

# point-scale cross-variogram
aod.list <- list(aod3k=aod3k.d, aod10=aod10.d)
vgm.ck <- deconvPointVgmForCoKriging(aod.list, model="Exp", ngroup=12, rd=0.75,
                                    fixed.range = 9e4)

## part 3: area-to-area Kriging prediction ----
ataStartCluster()
pred.ataok <- ataKriging(aod10.d, grid.pred, vgm.ck$aod10, showProgress = TRUE)
pred.ataok_combine <- ataKriging(aod.combine, grid.pred, vgm.ok_combine,
                                 showProgress = TRUE)
pred.atack <- ataCoKriging(aod.list, unknownVarId="aod3k", unknown=grid.pred,
                           ptVgms=vgm.ck, oneCondition=TRUE, auxRatioAdj=TRUE, showProgress = TRUE)

# convert result to raster
pred.ataok$pred <- exp(pred.ataok$pred)
pred.ataok_combine$pred <- exp(pred.ataok_combine$pred)
pred.atack$pred <- exp(pred.atack$pred)
pred.ataok.r <- rasterFromXYZ(pred.ataok[,2:4])
pred.ataok_combine.r <- rasterFromXYZ(pred.ataok_combine[,2:4])
pred.atack.r <- rasterFromXYZ(pred.atack[,2:4])

# display
pred <- stack(pred.ok.r, pred.ataok_combine.r, pred.ataok.r, pred.atack.r)
names(pred) <- c("point ok","ok_combine","ataok","atack")
spplot(pred)


## part 4: cross-validation ----
cv.ok <- c()
ext1 <- extent(aod3k)
ext2 <- extent(aod10)
xrange <- c(max(ext1[1],ext2[1]), min(ext1[2],ext2[2]))
for (i in 1:10) {
  print(i)
  x <- xrange[1] + (xrange[2]-xrange[1]) / 10 * c(i-1, i)
  cols <- colFromX(aod3k, x)
  aod3k.sub <- aod3k[,cols[1]:cols[2], drop=FALSE]

  cols <- colFromX(aod10, x)
  cols <- c(1:ncol(aod10))[-c(cols[1]:cols[2])]
  aod10.sub <- aod10[,cols,drop=FALSE]

  aod3k.sub.d <- discretizeRaster(aod3k.sub, 1500)
  aod10.sub.d <- discretizeRaster(aod10.sub, 1500)
  aod3k.sub.d$areaValues$value <- log(aod3k.sub.d$areaValues$value)
  aod10.sub.d$areaValues$value <- log(aod10.sub.d$areaValues$value)

  predi <- ataKriging(aod10.sub.d, aod3k.sub.d, vgm.ck$aod10, nmax=10)
  predi$value <- aod3k.sub.d$areaValues$value
  cv.ok <- rbind(cv.ok, predi)
}

cv.ok_combine <- c()
for (i in 1:10) {
  print(i)
  x <- xrange[1] + (xrange[2]-xrange[1]) / 10 * c(i-1, i)
  indx <- aod.combine$areaValues$centx >= x[1] & (aod.combine$areaValues$centx < x[2])
  areaId <-  aod.combine$areaValues$areaId[indx]

  unknown <- subsetDiscreteArea(aod.combine, areaId)
  predi <- ataKriging(subsetDiscreteArea(aod.combine, areaId, revSel = TRUE),
                      unknown, vgm.ok_combine, nmax = 10)
  cv.ok_combine <- rbind(cv.ok_combine, merge(unknown$areaValues, predi, by="areaId"))
}
cv.ok_combine <- cv.ok_combine[grepl("y_", cv.ok_combine$areaId),]

cv.ck <- ataCoKriging.cv(aod.list, unknownVarId = "aod3k", nfold = 10, ptVgms = vgm.ck,
                         showProgress = TRUE)
ataStopCluster()

# cross-validation for gstat ordinary Kriging
n <- floor(length(aod3k.pt) / 10)
cv.pointOk <- c()
for (i in 1:10) {
  nfrom <- (i-1) * n + 1
  nto <- nfrom + n - 1
  cur <- gstat::krige(value~1, aod3k.pt[-c(nfrom:nto),], newdata=aod3k.pt[nfrom:nto,],
               model=vgm.ok, nmax=10)
  cur$value <- aod3k.pt$value[nfrom:nto]
  if(is.null(cv.pointOk)) {
    cv.pointOk <- cur
  } else {
    cv.pointOk <- rbind(cv.pointOk, cur)
  }
}

evalAccuracy <- function(value, pred) {
  return(c(cor = cor(value, pred), mae = mean(abs(value - pred)), rmse = sqrt(mean((value - pred)^2))))
}

evalAccuracy(exp(cv.pointOk$value), exp(cv.pointOk$var1.pred))
evalAccuracy(exp(cv.ok$value), exp(cv.ok$pred))
evalAccuracy(exp(cv.ok_combine$value), exp(cv.ok_combine$pred))
evalAccuracy(exp(cv.ck$value), exp(cv.ck$pred))
