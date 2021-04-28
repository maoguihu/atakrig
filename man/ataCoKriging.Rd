\name{ataCoKriging}
\alias{ataCoKriging}
\alias{atpCoKriging}
\alias{ataCoKriging.cv}

\title{
Area-to-area, area-to-point coKriging prediciton, cross-validation.
}
\description{
Area-to-area, area-to-point coKriging prediciton, cross-validation.
}
\usage{
ataCoKriging(x, unknownVarId, unknown, ptVgms, nmax = 10, longlat = FALSE,
    oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE,
    showProgress = FALSE, nopar = FALSE, clarkAntiLog = FALSE)

atpCoKriging(x, unknownVarId, unknown0, ptVgms, nmax = 10, longlat = FALSE,
    oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE,
    showProgress = FALSE, nopar = FALSE)

ataCoKriging.cv(x, unknownVarId, nfold = 10, ptVgms, nmax = 10, longlat = FALSE,
    oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE,
    showProgress = FALSE, nopar = FALSE, clarkAntiLog = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
discretized areas of all variables, each is a discreteArea object.
}
  \item{unknownVarId}{
variable name (charaster) defined in x for prediction.
}
  \item{unknown}{
a discreted discreteArea object or data.frame[areaId,ptx,pty,weight] to be predicted.
}
  \item{unknown0}{
for points prediction or data.frame[ptx,pty] (one point per row) to be predicted.
}
  \item{nfold}{
number of fold for cross-validation. for leave-one-out cross-validation, nfold = nrow(x[[unknownVarId]]$areaValues).
}
  \item{ptVgms}{
point-scale direct and cross variograms, ataKrigVgm object.
}
  \item{nmax}{
max number of neighborhoods used for interpolation.
}
  \item{longlat}{
coordinates are longitude/latitude or not.
}
  \item{oneCondition}{
only one contrained condition for all points and all variables, \\sum_{i=1}^n\\lambda_i +\\sum_{j=1}^m\\beta_j =1, assuming expected means of variables known and constant with the study area.
}
  \item{meanVal}{
expected means of variables for oneCondition coKriging, data.frame(varId,value). If missing, simple mean values of areas from x will be used instead.
}
  \item{auxRatioAdj}{
for oneCondition Kriging, adjusting the auxiliary variable residue by a ratio between the primary variable mean and auxiliary variable mean.
}
  \item{showProgress}{
show progress bar for batch interpolation (multi destination areas).
}
  \item{nopar}{
disable parallel process in the function even if ataEnableCluster() has been called, mainly for internal use.
}
  \item{clarkAntiLog}{
for log-transformed input data, whether the estimated value should be adjusted(i.e. exponentiation).
}
}
%%\details{}
\value{
estimated value of destination area and its variance.
}
%%\note{%%  ~~further notes~~}
\references{
Clark, I., 1998. Geostatistical estimation and the lognormal distribution. Geocongress. Pretoria, RSA., [online] Available from: http://kriging.com/publications/Geocongress1998.pdf.
Goovaerts, P., 2008. Kriging and semivariogram deconvolution in the presence of irregular geographical units. Mathematical Geosciences 40 (1): 101-128.
Isaaks, E. H., Srivastava, R. M., 1989. An introduction to applied geostatistics. New York, Oxford University Press.
}
\seealso{
\link{deconvPointVgmForCoKriging}, \link{deconvPointCrossVgm}, \link{ataKriging}
}
\examples{
\donttest{
library(atakrig)
library(raster)

## demo data ----
rpath <- system.file("extdata", package="atakrig")
aod3k <- raster(file.path(rpath, "MOD04_3K_A2017042.tif"))
aod10 <- raster(file.path(rpath, "MOD04_L2_A2017042.tif"))

aod3k.d <- discretizeRaster(aod3k, 1500)
aod10.d <- discretizeRaster(aod10, 1500)
grid.pred <- discretizeRaster(aod3k, 1500, type = "all")

aod3k.d$areaValues$value <- log(aod3k.d$areaValues$value)
aod10.d$areaValues$value <- log(aod10.d$areaValues$value)

## area-to-area Kriging ----
# point-scale variogram from combined AOD-3k and AOD-10
aod.combine <- rbindDiscreteArea(aod3k.d, aod10.d)
sv.ok_combine <- deconvPointVgm(aod.combine, model="Exp", ngroup=12, rd=0.75)

# point-scale cross-variogram
aod.list <- list(aod3k=aod3k.d, aod10=aod10.d)
sv.ck <- deconvPointVgmForCoKriging(aod.list, model="Exp", ngroup=12, rd=0.75,
                                    fixed.range = 6.3e4)

# prediction
ataStartCluster(2) # parallel with 2 nodes
pred.ataok <- ataKriging(aod10.d, grid.pred, sv.ck$aod10, showProgress = TRUE)
pred.ataok_combine <- ataKriging(aod.combine, grid.pred, sv.ok_combine,
                                 showProgress = TRUE)
pred.atack <- ataCoKriging(aod.list, unknownVarId="aod10", unknown=grid.pred,
                           ptVgms=sv.ck, oneCondition=TRUE, auxRatioAdj=TRUE, showProgress = TRUE)
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
}
}
