\name{deconvPointVgm}
\alias{deconvPointVgm}
\alias{deconvPointCrossVariogram}
\alias{deconvPointVgmForCoKriging}
\title{
Point-scale variogram, cross-variogram deconvolution.
}
\description{
Point-scale variogram, cross-variogram deconvolution.
}
\usage{
deconvPointVgm(x, model = "Exp", maxIter = 100, fit.nugget = FALSE, fixed.range = NA, longlat = FALSE, maxSampleNum = 100, fig = TRUE, ...)
deconvPointCrossVariogram(x, y, xPointVgm, yPointVgm, model = "Exp", maxIter = 100, fit.nugget = FALSE, fixed.range = NA, longlat = FALSE, maxSampleNum = 100, fig = TRUE, ...)
deconvPointVgmForCoKriging(x, model = "Exp", maxIter = 100, fixed.range = NA, maxSampleNum = 100, fig = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
for \code{deconvPointVgm} and \code{deconvPointCrossVariogram}, x is a \code{discreteArea} object.

for \code{deconvPointVgmForCoKriging}, x is a list of \code{discreteArea} objects of all variables.

}
  \item{xPointVgm, yPointVgm}{
  point-scale variograms of x and y respectively, gstat variogramModel.
}
  \item{model}{
variogram model defined in \code{gstat::vgm()}, e.g. "Exp", "Sph", "Gau".
}
  \item{maxIter}{
max iteration number of deconvolution.
}
  \item{fit.nugget}{
fit variogram nugget or not.
}
  \item{fixed.range}{
variogram range fixed or not.
}
  \item{longlat}{
indicator whether coordinates are longitude/latitude.
}
  \item{maxSampleNum}{
to save memory and to reduce calculation time, for large number of discretized areas, only a number (maxSampleNum) of random sample will be used.
}
  \item{fig}{
whether to plot deconvoluted variogram.
}
  \item{\dots}{
additional paramters passed to \link{autofitVgm}.
}
}
\details{
The deconvolution algorithm is implemented according to Pierre Goovaerts, Math. Geosci., 2008, 40: 101-128.
}
\value{
  \item{pointVariogram}{
deconvoluted point variogram.
}
  \item{areaVariogram}{
fitted area variogram from area centroids.
}
  \item{experientialAreaVariogram}{
experiential area variogram from area centroids.
}
  \item{regularizedAreaVariogram}{
regularized area variogram from discretized area points and point variogram.
}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
