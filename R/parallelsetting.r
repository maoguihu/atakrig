## atakriging
## Author: Maogui Hu.

ataStartCluster <- function(spec = min(parallel::detectCores(), 8), ...) {
  cl <- getOption("ataKrigCluster")
  if(!is.null(cl)) try(snow::stopCluster(cl), silent = TRUE)

  cl <- snow::makeCluster(spec = spec, ...)
  doSNOW::registerDoSNOW(cl)
  options(ataKrigCluster = cl)
  cl
}


ataStopCluster <- function() {
  cl <- getOption("ataKrigCluster")
  if(!is.null(cl)) try(stopCluster(cl), silent = TRUE)
  options(ataKrigCluster = NULL)
}


ataClusterClearObj <- function() {
  cl <- getOption("ataKrigCluster")
  if(!is.null(cl)) try(clusterEvalQ(cl, "rm(list=ls())"), silent = TRUE)
}


ataIsClusterEnabled <- function() {
  return(!is.null(getOption("ataKrigCluster")))
}
