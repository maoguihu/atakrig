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
plotDeconvVgm(v, main = NULL, posx = NULL, posy = NULL, lwd = 2, showRegVgm = TRUE)
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{deconvPointVgmForCoKriging}, \link{deconvPointCrossVariogram}, \link{deconvPointCrossVariogram}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
