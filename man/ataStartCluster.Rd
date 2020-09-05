\name{ataStartCluster}
\alias{ataStartCluster}
\alias{ataStopCluster}
\alias{ataIsClusterEnabled}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Start/stop cluster parallel calculation.
}
\description{
Start/stop cluster parallel calculation for time consuming prediction.
ataIsClusterEnabled queries if cluster connections have been started by ataStartCluster.
}
\usage{
ataStartCluster(spec = min(parallel::detectCores(), 8), ...)
ataStopCluster()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
A specification appropriate to the type of cluster. See snow::\link[snow]{makeCluster}. By default, a maximum number of 8 slaves nodes can be creates on the local machine.
}
  \item{\dots}{
cluster type and option specifications.
}
}
