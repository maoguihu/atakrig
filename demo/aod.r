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
vgm.ok <- autofitVgm(aod3k.pt, ngroup=12, rd=0.75, fig = T)$model
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
vgm.ok_combine <- deconvPointVgm(aod.combine, model="Exp", ngroup=12, rd=0.75,
                                fit.nugget = TRUE)

# point-scale cross-variogram
aod.list <- list(aod3k=aod3k.d, aod10=aod10.d)
vgm.ck <- deconvPointVgmForCoKriging(aod.list, model="Exp", ngroup=12, rd=0.75,
                                    fixed.range = 9e4)

## part 3: area-to-area Kriging prediction ----
ataStartCluster()
pred.ataok <- ataKriging(aod10.d, grid.pred, vgm.ck$aod10, showProgress = T)
pred.ataok_combine <- ataKriging(aod.combine, grid.pred, vgm.ok_combine,
                                 showProgress = T)
pred.atack <- ataCoKriging(aod.list, unknownVarId="aod3k", unknown=grid.pred,
                           ptVgms=vgm.ck, oneCondition=T, auxRatioAdj=T, showProgress = T)

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
cv.ok <- ataKriging.cv(aod10.d, nfold = 10, ptVgm = vgm.ck$aod10, showProgress = T)
cv.ok_combine <- ataKriging.cv(aod.combine, nfold = 10, ptVgm = vgm.ok_combine, showProgress = T)
cv.ck <- ataCoKriging.cv(aod.list, unknownVarId = "aod3k", nfold = 10, ptVgms = vgm.ck,
                         showProgress = T)
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

# select pixels corresponding to original AOD3k
cv.ok_combineSub <- cv.ok_combine[(substr(cv.ok_combine$areaId, 1, 1) == "x"),]

evalAccuracy(exp(cv.pointOk$value), exp(cv.pointOk$var1.pred))
evalAccuracy(exp(cv.ok$value), exp(cv.ok$pred))
evalAccuracy(exp(cv.ck$value), exp(cv.ck$pred))
evalAccuracy(exp(cv.ok_combine$value), exp(cv.ok_combine$pred))
