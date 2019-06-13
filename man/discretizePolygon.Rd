\name{discretizePolygon}
\alias{discretizePolygon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discretize spatial polygons to points.
}
\description{
Discretize spatial polygons to points.
}
\usage{
discretizePolygon(x, cellsize, id = NULL, value = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a SpatialPolygonsDataFrame object.
}
  \item{cellsize}{
cell size of discretized grid.
}
  \item{id}{
unique polygon id.
}
  \item{value}{
polygon value.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
Point weight is normalized for each polygon. Weight need not to be the same for all points of a polygon. They can be assigned according to specific variables, such as population distribution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{discretizeRaster}
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
d1 <- discretizePolygon(observations, cellsize = 1500, id = "ID", value = "logObs")

plot(observations)
points(d1$discretePoints$ptx, d1$discretePoints$pty, pch=20)

# point-scale semivariogram deconvolution
gokDeconv <- deconvPointVgm(d1, model="Exp", ngroup=12, rd=0.75, fig=T)

# cross validation
cv <- ataKriging.cv(d1, nfold = nrow(d1$areaValues), gokDeconv$pointVariogram, nmax = 10)
cor(exp(cv$pred), exp(cv$value))

# prediction
d2 <- discretizePolygon(predictionLocations, cellsize = 1500, id = "ID")
pred <- ataKriging(d1, d2$discretePoints, gokDeconv$pointVariogram, nmax = 10)
pred$pred <- exp(pred$pred)
pred$var <- exp(pred$var)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
