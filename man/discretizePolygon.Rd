\name{discretizePolygon}
\alias{discretizePolygon}
\title{
Discretize spatial polygons to points.
}
\description{
Discretize spatial polygons to points.
}
\usage{
discretizePolygon(x, cellsize, id=NULL, value=NULL, showProgressBar=FALSE)
}
\arguments{
  \item{x}{
a SpatialPolygonsDataFrame object.
}
  \item{cellsize}{
cell size of discretized grid.
}
  \item{id}{
unique polygon id. if not given, polygons will be numbered from 1 to n accroding the record order.
}
  \item{value}{
polygon value. if not given, NA value will be assigned.
}
  \item{showProgressBar}{
whether show progress.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{discretizeRaster}, \link{ataKriging}
}
