\name{deconvPointVgm}
\alias{deconvPointVgm}
\alias{deconvPointCrossVgm}
\alias{deconvPointVgmForCoKriging}
\title{
Point-scale variogram, cross-variogram deconvolution.
}
\description{
Point-scale variogram, cross-variogram deconvolution.
}
\usage{
deconvPointVgm(x, model = "Exp", maxIter = 100,
    fixed.range = NA, longlat = FALSE, maxSampleNum = 100, fig = TRUE, ...)
deconvPointCrossVgm(x, y, xPointVgm, yPointVgm, model = "Exp",
    maxIter = 100, fixed.range = NA, longlat = FALSE,
    maxSampleNum = 100, fig = TRUE, ...)
deconvPointVgmForCoKriging(x, model = "Exp", maxIter = 100,
    fixed.range = NA, maxSampleNum = 100, fig = TRUE, ...)
}
\arguments{
  \item{x, y}{
for \code{deconvPointVgm} and \code{deconvPointCrossVgm}, x is a \code{discreteArea} object.

for \code{deconvPointVgmForCoKriging}, x is a list of \code{discreteArea} objects of all variables.

}
  \item{xPointVgm, yPointVgm}{
  point-scale variograms of x and y respectively, gstat variogramModel.
}
  \item{model}{
  commonly used variogram models supported, "Exp" for exponential model, "Sph" for spherical model, "Gau" for gaussian model.
%%variogram model defined in \code{gstat::vgm()}, e.g. "Exp", "Sph", "Gau".
}
  \item{maxIter}{
max iteration number of deconvolution.
}
  \item{fixed.range}{
a fixed variogram range for deconvoluted point-scale variogram.
}
  \item{longlat}{
indicator whether coordinates are longitude/latitude.
}
  \item{maxSampleNum}{
to save memory and to reduce calculation time, for large number of discretized areas, a number (maxSampleNum) of random sample will be used. The samples are collected by system sampling method.
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
Goovaerts, P., 2008. Kriging and semivariogram deconvolution in the presence of irregular geographical units. Mathematical Geosciences 40 (1): 101-128.
}
\seealso{
\link{ataKriging},\link{ataCoKriging}
}
\examples{
\donttest{
library(atakrig)
library(raster)

rpath <- system.file("extdata", package="atakrig")
aod3k <- raster(file.path(rpath, "MOD04_3K_A2017042.tif"))
aod3k <- crop(aod3k, extent(aod3k, 10,15,20,25))

aod3k.d <- discretizeRaster(aod3k, 1500)
grid.pred <- discretizeRaster(aod3k, 1500, type = "all")

sv.ok <- deconvPointVgm(aod3k.d, model="Exp", ngroup=12, rd=0.8, fig = FALSE)
#pred.ataok <- ataKriging(aod3k.d, grid.pred, sv.ok, showProgress = FALSE)
}
\donttest{
library(atakrig)
library(rgdal)

## load demo data from rtop package
#if (!require("rtop", quietly = TRUE)) message("rtop library is required for demo data.")
rpath <- system.file("extdata", package="rtop")
observations <- readOGR(rpath, "observations")

## point-scale variogram
obs.discrete <- discretizePolygon(observations, cellsize=1500, id="ID", value="obs")
pointsv <- deconvPointVgm(obs.discrete, model="Exp", ngroup=12, rd=0.75, fig=TRUE)
}
}
