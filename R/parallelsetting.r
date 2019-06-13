## atakriging
## Author: Maogui Hu, 2019.02.28

ataEnableCluster <- function(spec = min(parallel::detectCores(), 12), useDoSnow = TRUE, ...) {
  cl <- getOption("ataKrigCluster")
  if(!is.null(cl)) try(stopCluster(cl), silent = TRUE)

  if(useDoSnow) {
    if (!require(doSNOW)) {
      stop("failed to load doSNOW!")
    }
    cl <- makeCluster(spec = spec, ...)
    registerDoSNOW(cl)
  } else {
    if(!require(doParallel)) {
      stop("failed to load doParallel!")
    }
    cl <- makeCluster(spec = spec, ...)
    registerDoParallel(cl)
  }

  options(ataKrigCluster = cl)
  return(cl)
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
