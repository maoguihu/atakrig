\name{ataKriging}
\alias{ataKriging}
\alias{atpKriging}
\alias{ataKriging.cv}
\title{
Area-to-area, area-to-point ordinary kriging prediciton, cross-validation.
}
\description{
Area-to-area, area-to-point ordinary kriging prediciton, cross-validation.
}
\usage{
ataKriging(x, unknown, ptVgm, nmax = 10, longlat = FALSE, showProgress = TRUE, nopar = FALSE)
atpKriging(x, unknown0, ptVgm, nmax = 10, longlat=FALSE, showProgress=TRUE, nopar = FALSE)
ataKriging.cv(x, nfold = 10, ptVgm, nmax=10, longlat = FALSE, showProgress=TRUE, nopar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a discreteArea object: list(areaValues, discretePoints), where areaValues = data.frame(areaId,centx,centy,value); discretePoints = data.frame(areaId,ptx,pty,weight).
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
\link{deconvPointVgm}, \link{ataCoKriging}
}
\examples{
library(rgdal)

# load demo data from rtop package
rpath <- system.file("extdata", package="rtop")
observations <- readOGR(rpath, "observations")
predictionLocations <- readOGR(rpath, "predictionLocations")

hist(observations$obs)
observations$logObs <- log(observations$obs)

# discretize polygons
d1 <- discretizePolygon(observations, cellsize = 1500, idfield = "ID", valuefield = "logObs")

plot(observations)
points(d1$discretePoints$ptx, d1$discretePoints$pty, pch=20)

# point-scale semivariogram deconvolution
gokDeconv <- deconvPointVgm(d1, model="Exp", ngroup=12, rd=0.75, fig=T)

# cross validation
cv <- ataKriging.cv(d1, nfold = nrow(d1$areaValues), gokDeconv$pointVariogram, nmax = 10)
cor(exp(cv$pred), exp(cv$value))

# prediction
d2 <- discretizePolygon(predictionLocations, cellsize = 1500, idfield = "ID")
pred <- ataKriging(d1, d2$discretePoints, gokDeconv$pointVariogram, nmax = 10)
pred$pred <- exp(pred$pred)
pred$var <- exp(pred$var)

}
