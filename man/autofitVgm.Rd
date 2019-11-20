\name{autofitVgm}
\alias{autofitVgm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auto fit variogram for points.
}
\description{
Auto fit variogram for points.
}
\usage{
autofitVgm(x, y = x, ngroup = c(12, 15), rd = seq(0.3, 0.9, by = 0.1),
    model = c("Sph", "Exp", "Gau"), fit.nugget = TRUE, fixed.range = NA,
    longlat = FALSE, fig = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
values of areas, data.frame(areaId,centx,centy,value).
}
  \item{ngroup}{
number of bins to average from semivariogram cloud.
}
  \item{rd}{
ratio of max distance between points to be considered for bins.
}
  \item{model}{
variogram model defined in gstat::vgms(), e.g. "Exp", "Sph", "Gau".
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
  \item{fig}{
whether to plot fitted variogram.
}
  \item{\dots}{
additional parameters passed to gstat::vgm().
}
}
\value{
  \item{model}{fitted variogramModel.}
  \item{sserr}{fit error.}
  \item{bins}{binned gstatVariogram.}
}
\note{
The auto-search strategy was derived from \code{automap::autofitVariogram()}. The function tries different initial values of vgm to find the best fitted model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
