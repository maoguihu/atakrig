\name{ataEnableCluster}
\alias{ataEnableCluster}
\alias{ataStopCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enable/disable parallel calculation.
}
\description{
Enable/disable parallel calculation for some time consuming functions.
}
\usage{
ataEnableCluster(spec = max(parallel::detectCores(), 2), useDoSnow = TRUE, ...)
ataStopCluster()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
A specification appropriate to the type of cluster. See parallel::\link[parallel]{makeCluster}.
}
  \item{useDoSnow}{
use doSnow or doParallel for parallel calculation. Currently we recommend doSnow because a progressbar can be correctly displayed during calculation.
}
  \item{\dots}{
Passed to parallel::\link[parallel]{makeCluster} or snow::\link[snow]{makeCluster}.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
