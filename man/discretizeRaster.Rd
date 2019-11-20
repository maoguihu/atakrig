\name{discretizeRaster}
\alias{discretizeRaster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discretize raster to points.
}
\description{
Discretize raster to points.
}
\usage{
discretizeRaster(x, cellsize, type = "value", psf = "equal", sigma = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a RasterLayer object.
}
  \item{cellsize}{
cell size of discretized grid.
}
  \item{type}{
"value", "nodata", "all": whether only valid pixels, or only NODATA pixles, or all pixels extracted.
}
  \item{psf}{
PSF type, "equal", "gau", or user defined PSF matrix (normalized).
}
  \item{sigma}{
standard deviation for Gaussian PSF.
}
}
\value{
a discreteArea object: list(areaValues, discretePoints).
  \item{areaValues}{
values of areas: data.frame(areaId,centx,centy,value), where areaId is polygon id; centx, centy are centroids of polygons.
  }
  \item{discretePoints}{
discretized points of areas: data.frame(areaId,ptx,pty,weight), where ptx, pty are discretized points; by default, weight is equal for all points.
  }
}
\note{
Point weight is normalized for each polygon. Weight need not to be the same for all points of a polygon. They can be assigned according to specific variables, such as population distribution.
}
\seealso{
\link{discretizePolygon}, \link{ataCoKriging}
}
