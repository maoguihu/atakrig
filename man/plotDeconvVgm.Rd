\name{plotDeconvVgm}
\alias{plotDeconvVgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot deconvoluted point variogram.
}
\description{
Plot deconvoluted point variogram.
}
\usage{
plotDeconvVgm(v, main = NULL, posx = NULL, posy = NULL, lwd = 2, showRegVgm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
deconvoluted variogram, \code{ataKrigVgm}
}
  \item{main}{
title
}
  \item{posx, posy}{
position of legend
}
  \item{lwd}{
line width.
}
  \item{showRegVgm}{
show regularized area-scale variogram line or not.
}
}
\seealso{
\link{deconvPointVgmForCoKriging}, \link{deconvPointVgm}, \link{deconvPointCrossVgm}
}
