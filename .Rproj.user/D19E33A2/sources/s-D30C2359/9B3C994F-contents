library(atakrig)
library(raster)
library(gstat)

## demo data ----
rpath <- system.file("extdata", package="atakrig")
aod3k <- raster(file.path(rpath, "MOD04_3K_A2017042.tif"))
aod10 <- raster(file.path(rpath, "MOD04_L2_A2017042.tif"))

aod3k.d <- discretizeRaster(aod3k, 1000)
aod10.d <- discretizeRaster(aod10, 1000)
grid.pred <- discretizeRaster(aod3k, 1000, type = "all")

aod3k.d$areaValues$value <- log(aod3k.d$areaValues$value)
aod10.d$areaValues$value <- log(aod10.d$areaValues$value)

## gstat ok ----
aod3k.pt <- aod3k.d$areaValues
gok <- autofitVgm(aod3k.pt, ngroup=12, rd=0.75, fig = T)
unknown.pt <- grid.pred$areaValues
coordinates(aod3k.pt) = ~centx+centy
coordinates(unknown.pt) = ~centx+centy
pred.ok <- krige(value~1, aod3k.pt, newdata=unknown.pt, model=gok$model, nmax=10)

pred.ok$var1.pred <- exp(pred.ok$var1.pred)
pred.ok$var1.var <- exp(pred.ok$var1.var)
gridded(pred.ok) <- TRUE
spplot(pred.ok['var1.pred'])


## area-to-area kriging ----

set.seed(100)

# point-scale variogram from combined AOD-3k and AOD-10
aod.combine <- rbindDiscreteArea(aod3k.d, aod10.d)
sv.ok_combine <- deconvPointVgm(aod.combine, model="Exp", ngroup=12, rd=0.75, maxSampleNum=100)

# point-scale cross-variogram
aod.list <- list(aod3k=aod3k.d, aod10=aod10.d)
sv.ck <- deconvPointVgmForCoKriging(aod.list, model="Exp", ngroup=12, rd=0.75, fixed.range=8e4, maxSampleNum=100)

# prediction
ataEnableCluster()
pred.ataok <- ataKriging(aod10.d, grid.pred, sv.ck$aod10)
pred.ataok_combine <- ataKriging(aod.combine, grid.pred, sv.ok_combine)
pred.atack <- ataCoKriging(aod.list, unknownVar="aod10", unknown=grid.pred, ptVgms=sv.ck, oneCondition=TRUE, auxRatioAdj=T)
ataStopCluster()

# reverse log transform
pred.ataok$pred <- exp(pred.ataok$pred)
pred.ataok$var <- exp(pred.ataok$var)
pred.ataok_combine$pred <- exp(pred.ataok_combine$pred)
pred.ataok_combine$var <- exp(pred.ataok_combine$var)
pred.atack$pred <- exp(pred.atack$pred)
pred.atack$var <- exp(pred.atack$var)

# convert result to raster
pred.ataok.r <- rasterFromXYZ(pred.ataok[,-1])
pred.ataok_combine.r <- rasterFromXYZ(pred.ataok_combine[,-1])
pred.atack.r <- rasterFromXYZ(pred.atack[,-1])

# display
pred <- stack(aod3k, pred.ataok_combine.r$pred, pred.ataok.r$pred, pred.atack.r$pred)
names(pred) <- c("aod3k","ok_combine","ataok","atack")
spplot(pred)
