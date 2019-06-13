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
ataCoKriging(x, destVarId, unknown, ptVgms, nmax = 10, longlat = FALSE, oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE, showProgress = TRUE, nopar = FALSE)
atpCoKriging(x, destVarId, unknown0, ptVgms, nmax = 10, longlat = FALSE, oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE, showProgress = TRUE, nopar = FALSE)
ataCoKriging.cv(x, destVarId, nfold = 10, ptVgms, nmax = 10, longlat = FALSE, oneCondition = FALSE, meanVal = NULL, auxRatioAdj = TRUE, showProgress = TRUE, nopar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
discretized areas of all variables, list(
`varId1`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
`varId2`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
...,
`varIdn`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)))
}
  \item{destVarId}{
variable name (charaster) defined in x for prediction.
}
  \item{unknown}{
a discreted discreteArea object, or just data.frame(areaId,ptx,pty,weight).
}
  \item{unknown0}{
for points prediction, data.frame(ptx,pty), one point per row.
}
  \item{nfold}{
number of fold for cross-validation. for leave-one-out cross-validation, nfold = nrow(x[[destVarId]]$areaValues).
}
  \item{ptVgms}{
point-scale direct and cross variograms, ataKrigVgm.
}
  \item{nmax}{
max number of neighborhoods used for interpolation.
}
  \item{longlat}{
coordinates are longitude/latitude or not.
}
  \item{oneCondition}{
only one contrained condition for all points and all variables, ∑_(i=1)^nλ_i +∑_(j=1)^mβ_j =1, assuming expected means of variables known and constant with the study area.
}
  \item{meanVal}{
expected means of variables for oneCondition cokriging, data.frame(varId,value). If missing, simple mean values of areas from x will be used instead.
}
  \item{auxRatioAdj}{
for oneCondition kriging, adjusting the auxiliary variable residue by a ratio between the primary variable mean and auxiliary variable mean.
}
  \item{showProgress}{
show progress bar for batch interpolation (multi destination areas).
}
  \item{nopar}{
disable parallel process in the function even if ataEnableCluster() has been called, mainly for  internal use.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
estimated value of destination area and its variance.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{deconvPointVgmForCoKriging}, \link{deconvPointCrossVariogram}, \link{ataKriging}
}
\examples{
library(atakrig)
library(raster)

## demo data ----
rpath <- system.file("extdata", package="atakrig")
aod3k <- raster(file.path(rpath, "MOD04_3K_A2017042.tif"))
aod10 <- raster(file.path(rpath, "MOD04_L2_A2017042.tif"))

aod3k.d <- discretizeRaster(aod3k, 1000)
aod10.d <- discretizeRaster(aod10, 1000)
grid.pred <- discretizeRaster(aod3k, 1000, type = "all")

aod3k.d$areaValues$value <- log(aod3k.d$areaValues$value)
aod10.d$areaValues$value <- log(aod10.d$areaValues$value)

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
