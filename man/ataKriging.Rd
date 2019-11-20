\name{ataKriging}
\alias{ataKriging}
\alias{atpKriging}
\alias{ataKriging.cv}
\title{
Area-to-area, area-to-point ordinary Kriging prediciton, cross-validation.
}
\description{
Area-to-area, area-to-point ordinary Kriging prediciton, cross-validation.
}
\usage{
ataKriging(x, unknown, ptVgm, nmax = 10, longlat = FALSE,
    showProgress = TRUE, nopar = FALSE, clarkAntiLog = FALSE)
atpKriging(x, unknown0, ptVgm, nmax = 10, longlat=FALSE,
    showProgress=TRUE, nopar = FALSE)
ataKriging.cv(x, nfold = 10, ptVgm, nmax=10, longlat = FALSE,
    showProgress=TRUE, nopar = FALSE, clarkAntiLog = FALSE)
}
\arguments{
  \item{x}{
a discreteArea object: list(areaValues, discretePoints), where
    areaValues: data.frame(areaId,centx,centy,value)
    discretePoints: data.frame(areaId,ptx,pty,weight)
}
  \item{unknown}{
a discreted discreteArea object, or just data.frame(areaId,ptx,pty,weight).
}
  \item{unknown0}{
for points prediction, data.frame(ptx,pty), one point per row.
}
  \item{nfold}{
number of fold for cross-validation. for leave-one-out cross-validation, nfold = nrow(x$areaValues).
}
  \item{ptVgm}{
point scale variogram, ataKrigVgm.
}
  \item{nmax}{
max number of neighborhoods used for interpolation.
}
  \item{longlat}{
coordinates are longitude/latitude or not.
}
  \item{showProgress}{
show progress bar for batch interpolation (multi destination areas).
}
  \item{nopar}{
disable parallel process in the function even if \link{ataStartCluster}() has been called, mainly for  internal use.
}
  \item{clarkAntiLog}{
for log-transformed input data, whether the estimated value should be adjusted(i.e. exponentiation).
}
}
%%\details{}
\value{
estimated value of destination area and its variance.
}
\references{
Clark, I., 1998. Geostatistical estimation and the lognormal distribution. Geocongress. Pretoria, RSA., [online] Available from: http://kriging.com/publications/Geocongress1998.pdf.
Goovaerts, P., 2008. Kriging and semivariogram deconvolution in the presence of irregular geographical units. Mathematical Geosciences 40 (1): 101-128.
Isaaks, E. H., Srivastava, R. M., 1989. An introduction to applied geostatistics. New York, Oxford University Press.
}
\seealso{
\link{deconvPointVgm}, \link{ataCoKriging}
}
\examples{
\dontrun{
library(rgdal)

## load demo data from rtop package
if (!require("rtop", quietly = TRUE)) install.packages("rtop")
rpath <- system.file("extdata", package="rtop")
observations <- readOGR(rpath, "observations")

# hist(observations$obs)
observations$obs <- log(observations$obs)

## point-scale variogram
obs.discrete <- discretizePolygon(observations, cellsize=1500, id="ID", value="obs")
pointsv <- deconvPointVgm(obs.discrete, model="Exp", ngroup=12, rd=0.75, fig=TRUE)

## prediction
predictionLocations <- readOGR(rpath, "predictionLocations")
pred.discrete <- discretizePolygon(predictionLocations, cellsize = 1500, id = "ID")
pred <- ataKriging(obs.discrete, pred.discrete, pointsv$pointVariogram, nmax = 10)
pred$pred <- exp(pred$pred)
pred$var <- exp(pred$var)

## cross validation
pred.cv <- ataKriging.cv(obs.discrete, nfold=nrow(obs.discrete$areaValues), pointsv,
    showProgress = TRUE)
pred.cv$pred <- exp(pred.cv$pred)
pred.cv$obs <- exp(pred.cv$value)

# pearson correlation
cor(pred.cv$obs, pred.cv$pred)

# mae
mean(abs(pred.cv$obs - pred.cv$pred))

# rmse
sqrt(mean((pred.cv$obs - pred.cv$pred)^2))
}
}
