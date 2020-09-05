\name{ataSetNumberOfThreadsForOMP}
\alias{ataSetNumberOfThreadsForOMP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set number of threads for OpenMP.
}
\description{
Set number of threads for OpenMP.
}
\usage{
ataSetNumberOfThreadsForOMP(num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num}{
  An integer number of threads for OpenMP.
}
}
\details{
  The deconvolution of variogram is computation intensive. Some parts of them is coded by Rcpp with OpenMP enabled. By default, the number of threads created by OpenMP is the number of local machine cores.
  It should be noted that OpenMP is not supported for macOS since R 4.0.0.
}
\seealso{
\link{ataStartCluster}
}
