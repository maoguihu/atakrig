## atakrig
## Function: Deconvolution of point-scale variogram/cross-variogram.
## Author: Maogui Hu, 2019.02.28

# require(sp)
# require(gstat)

## deconvPointVgmForCoKriging ----
# Input:
#   x: discretized areas, list(
#      `varId1`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#      `varId2`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)),
#       ...,
#      `varIdn`=list(areaValues=data.frame(areaId,centx,centy,value), discretePoints=data.frame(areaId,ptx,pty,weight)))
#  model: variogram model type, e.g. "Exp", "Sph", "Gau", "Mat". Calling vgm() returns available models.
#  fig: plot deconvoluted variograms when finished.
# Output:
#   a list of direct/cross variograms.
deconvPointVgmForCoKriging <- function(x, model="Exp", maxIter=100, fixed.range=NA, maxSampleNum=100, fig=TRUE, ...) {
  varnames <- sort(names(x))
  vgms <- list()

  # vargioram
  for (i in 1:length(varnames)) {
    id <- varnames[i]
    if(!hasName(x[[id]], "discretePoints")) {
      x[[id]]$discretePoints <- cbind(x[[id]]$areaValues[,1:3], data.frame(weight=rep(1,nrow(x[[id]]$areaValues))))
      names(x[[id]]$discretePoints)[2:3] <- c("ptx","pty")
    }

    if(nrow(x[[id]]$areaValues) > maxSampleNum) {
      warnings("Too many points for generating point-pairs. Only a sample will be used.")
      indx <- sort(sample(nrow(x[[id]]$areaValues), maxSampleNum))
      x[[id]] <- subsetDiscreteArea(x[[id]], x[[id]]$areaValues[indx,1])
    }

    cat(sprintf("Deconvoluting variogram of %s ...\n", id))
    vgms[[id]] <- deconvPointVgm(x[[id]], model = model, maxIter = maxIter, fixed.range=fixed.range, fig=FALSE, ...)
    stopifnot(!is.null(vgms[[id]]))
    if(vgms[[id]]$status == 0) {
      warning(sprintf("deconvolution failed for %s!", id))
    }
  }

  # cross-variogram
  for (i in 1:(length(varnames)-1)) {
    id1 <- varnames[i]
    for (j in (i+1):length(varnames)) {
      id2 <- varnames[j]
      cat(sprintf("Deconvoluting cross-variogram between %s and %s ...\n", id1, id2))
      id <- .crossName(id1,id2)
      vgms[[id]] <- deconvPointCrossVariogram(x[[id1]], x[[id2]],
                                                   vgms[[id1]]$pointVariogram,
                                                   vgms[[id2]]$pointVariogram,
                                                   model = model, maxIter = maxIter, fixed.range=fixed.range,
                                                   fig=FALSE, ...)
      stopifnot(!is.null(vgms[[id]]))
      if(vgms[[id]]$status == 0) {
        warning(sprintf("deconvolution failed for %s!", id))
      }
    }
  }

  if (!is.na(fixed.range)) {
    posdef = function(X) {
      q = eigen(X)
      d = q$values
      d[d < 0] = 0
      q$vectors %*% diag(d, nrow = length(d)) %*% t(q$vectors)
    }

    psill = matrix(NA, nrow = length(varnames), ncol = length(varnames))
    for (i in 1:length(varnames)) {
      for (j in i:length(varnames)) {
        id = ifelse(i == j, varnames[i], .crossName(varnames[i], varnames[j]))
        psill[i, j] = psill[j, i] = vgms[[id]]$pointVariogram[1,"psill"]
      }
    }
    psill = posdef(psill)
    for (i in 1:length(varnames)) {
      for (j in i:length(varnames)) {
        id = ifelse(i == j, varnames[i], .crossName(varnames[i], varnames[j]))
        vgms[[id]]$pointVariogram[1,"psill"] = psill[i, j]
      }
    }
  }

  if(fig) plotDeconvVgm(vgms, main = "Deconvoluted variograms/cross-variograms")

  class(vgms) <- c("list", "ataKrigVgm")
  return(vgms)
}


## deconvPointVgm: Point scale variogram deconvolution according to Pierre Goovaerts, Math. Geosci., 2008, 40: 101-128. ----
# Input:
#   x: discretized area, list(areaValues,discretePoints):
#      areaValues: values of areas (data.frame: areaId,centx,centy,value).
#      discretePoints: discretized points of areas (data.frame: areaId,ptx,pty,weight); the weight is normalized.
#   model: variogram model type, e.g. "Exp", "Sph", "Gau", "Mat". Calling vgm() returns available models.
#   maxIter: max iteration number of deconvolution.
#   fit.nugget: fit variogram nugget or not.
#   fixed.range: variogram range fixed or not.
#	  longlat: indicator whether coordinates are longitude/latitude
#   maxSampleNum: to save memory and to reduce calculation time, for large number of discretized areas,
#      only a number (maxSampleNum) of random sample will be used.
#   fig: whether to plot deconvoluted variogram.
# Output:
#   list(pointVariogram,             # deconvoluted point variogram.
#        areaVariogram,              # fitted area variogram from area centroids.
#        experientialAreaVariogram,  # experiential area variogram from area centroids.
#        regularizedAreaVariogram    # regularized area variogram from discretized area points and point variogram.
#        ).
deconvPointVgm <- function(x, model="Exp", maxIter=100, fit.nugget=FALSE, fixed.range=NA, longlat=FALSE, maxSampleNum=100, fig=TRUE, ...) {
  # Pierre Goovaerts suggested a different weight (similar to fit.method = 2) in fitting the experimental variogram.

  if(!hasName(x, "discretePoints")) x$discretePoints <- cbind(x$areaValues[,1:3], data.frame(weight=rep(1,nrow(x$areaValues))))

  if(nrow(x$areaValues) > maxSampleNum) {
    warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sample(nrow(x$areaValues), maxSampleNum))
    x <- subsetDiscreteArea(x, x$areaValues[indx,1])
  }

  ## 1. area scale semivariogram
  areaSvmodel <- autofitVgm(x$areaValues, fit.nugget=fit.nugget, fixed.range=fixed.range, longlat=longlat, model=model, ...)
  if(is.null(areaSvmodel)) {
    cat("Fitting area-scale variogram model failed!\n")
    return(NULL)
  }
  if(nrow(x$areaValues) == nrow(x$discretePoints)) {
    vgms <- list(pointVariogram = areaSvmodel$model,
                  areaVariogram = areaSvmodel$model,
                  experientialAreaVariogram = areaSvmodel$bins,
                  regularizedAreaVariogram = NULL,
                  status = 1,
                  sserr = areaSvmodel$sserr)
    class(vgms) <- c("list", "ataKrigVgm")
    return(vgms)
  }
  ngroup <- areaSvmodel$ngroup
  rd <- areaSvmodel$rd

  ## 0. prepare
  x$areaValues <- x$areaValues[order(x$areaValues$areaId),]

  # distance between area centroids
  areaDistByCentroid <- spDists(as.matrix(x$areaValues[,2:3]), longlat=longlat)

  # entry index for each area
  uId <- sort(unique(x$discretePoints[,1]))

  # distances between discretized area points
  hasCluster <- !is.null(getOption("ataKrigCluster"))
  if(!hasCluster) {
    areaDistByPts <- list()
    for(i in 1:length(uId)) {
      areaDistByPts[[i]] <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId[i])
      for(j in i:length(uId)) {
        areaIndexJ <- which(x$discretePoints[,1] == uId[j])
        areaDistByPts[[i]][[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                           as.matrix(x$discretePoints[areaIndexJ,2:3]),
                                           longlat=longlat)
      }
    }
  } else {
    areaDistByPts <- foreach(i = 1:length(uId), .packages = "sp") %dopar%  {
      areaDistByPts <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId[i])
      for(j in i:length(uId)) {
        areaIndexJ <- which(x$discretePoints[,1] == uId[j])
        areaDistByPts[[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                           as.matrix(x$discretePoints[areaIndexJ,2:3]),
                                           longlat=longlat)
      }
      return(areaDistByPts)
    }
    ataClusterClearObj()
  }

  ## 2. initialize
  pointSvmodel <- areaSvmodel$model
  dgGroup <- areaSvmodel$bins
  gamaExp <- variogramLine(areaSvmodel$model, covariance=FALSE, dist_vector=dgGroup$dist)$gamma
  s2Exp <- areaSvmodel$model$psill

  ## 3. regularize
  gamaReg <- .calcSvDgGroup(.calcSvCloudAreaByPointVgm(x$discretePoints, pointSvmodel,
                                         areaDistByCentroid = areaDistByCentroid,
                                         areaDistByPts = areaDistByPts),
                           ngroup = ngroup, rd = rd)$gamma

  ## 4. difference of areal regularized and experimental semivariogram
  D0 <- mean(abs(gamaReg-gamaExp)/gamaExp)

  ## 5. backup for comparison
  pointSvmodelOpt <- pointSvmodel
  gamaOpt <- gamaReg
  DOpt <- D0

  nvib <- 0
  bNewW <- FALSE
  status <- -1
  for(iter in 1:maxIter) {
    cat(sprintf("\riterating: %d", iter))

    ## 6. update regularized semivariance
    if(!bNewW)
      wl <- 1 + (gamaExp-gamaOpt)/(s2Exp*sqrt(iter))
    gPoint <- variogramLine(pointSvmodelOpt, covariance=FALSE, dist_vector=dgGroup$dist)$gamma * wl

    ## 7. fit new point scale semivariogram
    dgGroup$gamma <- gPoint
    pointSvmodel <- .fitPointVgm(dgGroup, pointSvmodel[nrow(pointSvmodel),1], fit.nugget=fit.nugget, fixed.range=fixed.range)$model
    if(is.null(pointSvmodel)) {
      pointSvmodel <- pointSvmodelOpt
      status <- 0
      # warning("fitting point-scale variogram failed!")
      break
    }

    ## 8. regularize
    gamaReg <- .calcSvDgGroup(.calcSvCloudAreaByPointVgm(x$discretePoints, pointSvmodel,
                                           areaDistByCentroid = areaDistByCentroid,
                                           areaDistByPts = areaDistByPts),
                             ngroup = ngroup, rd = rd)$gamma

    ## 9. difference of areal regularized and experimental semivariogram
    D1 <- mean(abs(gamaReg-gamaExp)/gamaExp)

    # 10. stop criterion
    if(max(abs(wl-1)) < 0.001) {
      status <- 1
      break
    }
    if(abs(D1-DOpt)/DOpt <= 0.01) {
      nvib <- nvib+1
      if(nvib >= 5) {
        status <- 2
        break
      }
    } else {
      nvib <- 0
    }
    if(D1/D0 <= 0.01) {
      status <- 3
      break
    }

    if(D1 < DOpt) {
      pointSvmodelOpt <- pointSvmodel
      gamaOpt <- gamaReg
      DOpt <- D1
      bNewW <- FALSE
    } else {
      pointSvmodel <- pointSvmodelOpt
      wl <- 1 + (wl-1)/2
      bNewW <- TRUE
    }
  }
  if(iter == maxIter) {
    status <- 4
  }
  cat("\r",rep(" ",15),"\r")

  vgms <- list(pointVariogram = pointSvmodel,
                 areaVariogram = areaSvmodel$model,
                 experientialAreaVariogram = areaSvmodel$bins,
                 regularizedAreaVariogram = data.frame(dist=areaSvmodel$bins$dist, gamma=gamaOpt),
                 status = status,
                 sserr = areaSvmodel$sserr)
  class(vgms) <- c("list", "ataKrigVgm")

  if(fig) try(plotDeconvVgm(vgms), silent = T)

  return(vgms)
}


## deconvPointCrossVariogram: Point scale cross-variogram deconvolution. ----
# Input:
#   x, y: discretized areas, list(discretePoints, areaValues):
#       areaValues: values of areas, data.frame(areaId,centx,centy,value).
#       discretePoints: discretized points of areas, data.frame(areaId,ptx,pty,weight), the weight is normalized.
#   xPointVgm, yPointVgm: point-scale variograms of x and y respectively, gstat variogramModel.
#   model: variogram model type, e.g. "Exp", "Sph", "Gau", "Mat". Calling vgm() returns available models.
#   maxIter: max iteration number of deconvolution.
#	  longlat: indicator whether coordinates are longitude/latitude.
#   maxSampleNum: to save memory and to reduce calculation time, for large number of discretized areas,
#      only a number (maxSampleNum) of random sample will be used.
#   fig: whether to plot deconvoluted variogram.
# Output:
#   list(pointVariogram,             # deconvoluted point variogram.
#        areaVariogram,              # fitted area variogram from area centroids.
#        experientialAreaVariogram,  # experiential area variogram from area centroids.
#        regularizedAreaVariogram    # regularized area variogram from discretized area points and point variogram.
#        ).
deconvPointCrossVariogram <- function(x, y, xPointVgm, yPointVgm, model="Exp", maxIter=100, fit.nugget=FALSE, fixed.range=NA,
                                           longlat=FALSE, maxSampleNum=100, fig=TRUE, ...) {

  if(!hasName(x, "discretePoints")) x$discretePoints <- cbind(x$areaValues[,1:3], data.frame(weight=rep(1,nrow(x$areaValues))))
  if(!hasName(y, "discretePoints")) y$discretePoints <- cbind(y$areaValues[,1:3], data.frame(weight=rep(1,nrow(y$areaValues))))

  if(nrow(x$areaValues) > maxSampleNum) {
    warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sample(nrow(x$areaValues), maxSampleNum))
    x <- subsetDiscreteArea(x, x$areaValues[indx,1])
  }
  if(nrow(y$areaValues) > maxSampleNum) {
    warnings("Too many points for generating point-pairs. Only a sample will be used.")
    indx <- sort(sample(nrow(y$areaValues), maxSampleNum))
    y <- subsetDiscreteArea(y, y$areaValues[indx,1])
  }

  ## 1. area scale cross-semivariogram
  areaCrossVgm <- autofitVgm(x$areaValues, y$areaValues, fit.nugget=fit.nugget, fixed.range=fixed.range, longlat=longlat, model=model, ...)
  if(is.null(areaCrossVgm$model)) {
    cat("Fitting area-scale cross-variogram model failed!\n")
    return(NULL)
  }
  if(nrow(x$areaValues) == nrow(x$discretePoints) && nrow(y$areaValues) == nrow(y$discretePoints)) {
    vgms <- list(pointVariogram = areaCrossVgm$model,
                  areaVariogram = areaCrossVgm$model,
                  experientialAreaVariogram = areaCrossVgm$bins,
                  regularizedAreaVariogram = NULL,
                  status = 1,
                  sserr = areaCrossVgm$sserr)
    class(vgms) <- c("list", "ataKrigVgm")
    return(vgms)
  }
  ngroup <- areaCrossVgm$ngroup
  rd <- areaCrossVgm$rd

  ## 0. prepare
  x$areaValues <- x$areaValues[order(x$areaValues$areaId),]
  y$areaValues <- y$areaValues[order(y$areaValues$areaId),]

  # distance between area centroids
  areaDistByCentroid12 <- spDists(as.matrix(x$areaValues[,2:3]), as.matrix(y$areaValues[,2:3]), longlat=longlat)

  # entry index for each area
  uId1 <- sort(unique(x$discretePoints[,1]))
  uId2 <- sort(unique(y$discretePoints[,1]))

  # distances between discretized area points
  hasCluster <- !is.null(getOption("ataKrigCluster"))
  if(!hasCluster) {
    areaDistByPts1 <- list()
    for(i in 1:length(uId1)) {
      areaDistByPts1[[i]] <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in i:length(uId1)) {
        areaIndexJ <- which(x$discretePoints[,1] == uId1[j])
        areaDistByPts1[[i]][[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                            as.matrix(x$discretePoints[areaIndexJ,2:3]),
                                            longlat=longlat)
      }
    }
    areaDistByPts2 <- list()
    for(i in 1:length(uId2)) {
      areaDistByPts2[[i]] <- list()
      areaIndexI <- which(y$discretePoints[,1] == uId2[i])
      for(j in i:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPts2[[i]][[j]] <- spDists(as.matrix(y$discretePoints[areaIndexI,2:3]),
                                            as.matrix(y$discretePoints[areaIndexJ,2:3]),
                                            longlat=longlat)
      }
    }
    areaDistByPts12 <- list()
    for(i in 1:length(uId1)) {
      areaDistByPts12[[i]] <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in 1:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPts12[[i]][[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                             as.matrix(y$discretePoints[areaIndexJ,2:3]),
                                             longlat=longlat)
      }
    }
  } else {
    areaDistByPts1 <- foreach(i = 1:length(uId1), .packages = "sp") %dopar% {
      areaDistByPts1 <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in i:length(uId1)) {
        areaIndexJ <- which(x$discretePoints[,1] == uId1[j])
        areaDistByPts1[[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                            as.matrix(x$discretePoints[areaIndexJ,2:3]),
                                            longlat=longlat)
      }
      return(areaDistByPts1)
    }
    areaDistByPts2 <- foreach(i = 1:length(uId2), .packages = "sp") %dopar% {
      areaDistByPts2 <- list()
      areaIndexI <- which(y$discretePoints[,1] == uId2[i])
      for(j in i:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPts2[[j]] <- spDists(as.matrix(y$discretePoints[areaIndexI,2:3]),
                                            as.matrix(y$discretePoints[areaIndexJ,2:3]),
                                            longlat=longlat)
      }
      return(areaDistByPts2)
    }
    areaDistByPts12 <- foreach(i = 1:length(uId1), .packages = "sp") %dopar% {
      areaDistByPts12 <- list()
      areaIndexI <- which(x$discretePoints[,1] == uId1[i])
      for(j in 1:length(uId2)) {
        areaIndexJ <- which(y$discretePoints[,1] == uId2[j])
        areaDistByPts12[[j]] <- spDists(as.matrix(x$discretePoints[areaIndexI,2:3]),
                                             as.matrix(y$discretePoints[areaIndexJ,2:3]),
                                             longlat=longlat)
      }
      return(areaDistByPts12)
    }
    ataClusterClearObj()
  }

  ## 2. initialize
  pointCrossVgm <- areaCrossVgm$model
  dgGroup <- areaCrossVgm$bins
  gamaExp <- variogramLine(areaCrossVgm$model, covariance=FALSE, dist_vector=dgGroup$dist)$gamma
  s2Exp <- areaCrossVgm$model$psill

  ## 3. regularize
  crossDgScatter <- .calcCrossSvCloudAreaByPointVgm(x$discretePoints, y$discretePoints, xPointVgm, yPointVgm,
                                                   xyPointCrossVgm = pointCrossVgm, areaDistByCentroid12,
                                                   areaDistByPts1, areaDistByPts2, areaDistByPts12)
  gamaReg <- .calcSvDgGroup(crossDgScatter, ngroup = ngroup, rd = rd)$gamma

  ## 4. difference of areal regularized and experimental semivariogram
  D0 <- mean(abs(gamaReg-gamaExp)/gamaExp)

  ## 5. backup for comparison
  pointCrossVgmOpt <- pointCrossVgm
  gamaOpt <- gamaReg
  DOpt <- D0

  nvib <- 0
  bNewW <- FALSE
  status <- -1
  for(iter in 1:maxIter) {
    cat(sprintf("\riterating: %d", iter))

    ## 6. update regularized semivariance
    if(!bNewW)
      wl <- 1 + (gamaExp-gamaOpt)/(s2Exp*sqrt(iter))
    gPoint <- variogramLine(pointCrossVgmOpt, covariance=FALSE, dist_vector=dgGroup$dist)$gamma * wl

    ## 7. fit new point scale semivariogram
    dgGroup$gamma <- gPoint
    pointCrossVgm <- .fitPointVgm(dgGroup, pointCrossVgm[nrow(pointCrossVgm),1], fit.nugget=fit.nugget, fixed.range=fixed.range, ...)$model
    if(is.null(pointCrossVgm)) {
      pointCrossVgm <- pointCrossVgmOpt
      status <- 0
      # warning("fitting point-scale cross-variogram failed!")
      break
    }

    ## 8. regularize
    crossDgScatter <- .calcCrossSvCloudAreaByPointVgm(x$discretePoints, y$discretePoints, xPointVgm, yPointVgm,
                                                     xyPointCrossVgm = pointCrossVgm, areaDistByCentroid12,
                                                     areaDistByPts1, areaDistByPts2, areaDistByPts12)
    gamaReg <- .calcSvDgGroup(crossDgScatter, ngroup = ngroup, rd = rd)$gamma

    ## 9. difference of areal regularized and experimental semivariogram
    D1 <- mean(abs(gamaReg-gamaExp)/gamaExp)

    # 10. stop criterion
    if(max(abs(wl-1)) < 0.001) {
      status <- 1
      break
    }
    if(abs(D1-DOpt)/DOpt <= 0.01) {
      nvib <- nvib+1
      if(nvib >= 5) {
        status <- 2
        break
      }
    } else {
      nvib <- 0
    }
    if(D1/D0 <= 0.01) {
      status <- 3
      break
    }

    if(D1 < DOpt) {
      pointCrossVgmOpt <- pointCrossVgm
      gamaOpt <- gamaReg
      DOpt <- D1
      bNewW <- FALSE
    } else {
      pointCrossVgm <- pointCrossVgmOpt
      wl <- 1 + (wl-1)/2
      bNewW <- TRUE
    }
  }
  if(iter == maxIter) {
    status <- 4
  }
  cat("\r",rep(" ",15),"\r")

  vgms <- list(pointVariogram = pointCrossVgm,
                 areaVariogram = areaCrossVgm$model,
                 experientialAreaVariogram = areaCrossVgm$bins,
                 regularizedAreaVariogram = data.frame(dist=areaCrossVgm$bins$dist, gamma=gamaOpt),
                 status = status,
                 sserr = areaCrossVgm$sserr)
  class(vgms) <- c("list", "ataKrigVgm")

  if(fig) try(plotDeconvVgm(vgms), T)
  return(vgms)
}


## plotDeconvVgm: Plot deconvoluted point variogram. ----
plotDeconvVgm <- function(v, main=NULL, posx=NULL, posy=NULL, lwd=2, showRegVgm=TRUE) {
  .plotDeconvVgm <- function(v, main) {
    xlim <- c(0, 1.1 * max(v$experientialAreaVariogram$dist))
    xx <- seq(0, xlim[2], length=100)
    yy <- variogramLine(v$areaVariogram, covariance=FALSE, dist_vector=xx)$gamma
    yy2 <- variogramLine(v$pointVariogram, covariance=FALSE, dist_vector=xx)$gamma

    ylim <- c(min(0, min(v$experientialAreaVariogram$gamma)),
              1.1 * max(c(v$experientialAreaVariogram$gamma, yy, yy2)))
    plot(v$experientialAreaVariogram$dist, v$experientialAreaVariogram$gamma,
         xaxs="i", yaxs="i", col="blue", xlim=xlim, ylim=ylim, ann = F,
         xlab="distance", ylab="semivariance", main=main)
    lines(xx, yy, col="blue", lwd=lwd)
    lines(xx, yy2, col="red", lwd=lwd)
    if(showRegVgm)
      lines(v$regularizedAreaVariogram$dist, v$regularizedAreaVariogram$gamma, col="black", type="l", lty=2, lwd=lwd)
    mtext(main, side = 3, line = 0.2, cex = 0.8)
  }

  def.par <- par(no.readonly = TRUE)
  par(mar=c(2.5, 3.5, 1, .5)) # reduce the margins around the figure
  par(mgp=c(1.5, .5, 0)) #  reduce the spacing between the figure plotting region and the axis labels
  par(oma=c(1, 1, 2.5, 0)) #  add an outer margin to the top of the graph

  if (hasName(v, "regularizedAreaVariogram")) {
    .plotDeconvVgm(v, main = "")
    if(is.null(posx)) posx <- "bottomright"
  } else {
    m0 <- (sqrt(1+8*length(v))-1)/2
    n <- sort(names(v)[1:m0])
    m1 <- m0*(m0-1)/2 + m0
    m <- matrix(0, nrow=m0, ncol=m0)
    m[lower.tri(m, diag = T)] <- seq(m1)
    m[1,2:m0] <- m1 + 1
    layout(mat=m)

    for (i in 1:m0) {
      .plotDeconvVgm(v[[n[i]]], main = n[i])
      for (j in 1:m0) {
        if(j > i) {
          n2 <- .crossName(n[i],n[j])
          .plotDeconvVgm(v[[n2]], main = n2)
        }
      }
    }
    plot.new()
    if(is.null(posx)) posx <- "left"
  }

  if(showRegVgm) {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram", "Regularized area variogram"),
           pch=c(1,NA,NA,NA), lty=c(NA,1,1,2), lwd=c(NA,lwd,lwd,lwd),
           col=c("blue","blue","red","black"))
  } else {
    legend(posx, posy, bty="n",
           legend=c("Empirical area variogram", "Fitted area variogram",
                    "Deconvoluted point variogram"),
           pch=c(1,NA,NA), lty=c(NA,1,1), lwd=c(NA,lwd,lwd),
           col=c("blue","blue","red","black"))
  }

  if(is.null(main)) main <- "Deconvoluted variogram"
  mtext(main, outer = TRUE,  side = 3, cex = 1.2, line = 0)
  mtext("distance", outer = TRUE,  side = 1, cex = 1.05, line = 0)
  mtext("semivariance", outer = TRUE,  side = 2, cex = 1.05, line = -1)

  par(def.par)
}


## autofitVgm: Auto fit variogram for points or areas according to the centroids. ----
# Input:
#   x, y: values of areas, data.frame(areaId,centx,centy,value).
#   model: as.character(vgm()$short)[-1]
autofitVgm <- function(x, y=x, ngroup=c(12,15), rd=seq(0.3,0.9,by=0.1), model=c("Sph","Exp","Gau"),
                       fit.nugget=TRUE, fixed.range=NA, longlat=FALSE, fig=FALSE, ...) {
  areaSvCloud <- .calcSvCloud(x[,2:4], y[,2:4], longlat=longlat)
  mSel <- NULL
  for (i in 1:length(ngroup)) {
    for(j in 1:length(rd)) {
      areaVgm <- .calcSvDgGroup(areaSvCloud, ngroup=ngroup[i], rd=rd[j])
      if(nrow(areaVgm) < 10) next

      m <- .fitPointVgm(vgmdef=areaVgm, model=model, fit.nugget=fit.nugget, fixed.range=fixed.range, ...)
      if(is.null(m$model)) next

      m$ngroup <- ngroup[i]
      m$rd <- rd[j]
      if(is.null(mSel)) {
        mSel <- m
      } else {
        if(m$sserr < mSel$sserr) mSel <- m
      }
    }
  }

  if(fig && !is.null(mSel)) print(plot(mSel$bins, mSel$model))
  return(mSel)
}


## .fitPointVgm: [internal use only] Fit variogram for points. ----
#      vgmdef generted by .calcSvDgGroup()
.fitPointVgm <- function(vgmdef, model=c("Sph","Exp","Gau"), fit.nugget=TRUE, fixed.range=NA, fig=FALSE, ...) {
  init_nugget <- ifelse(!fit.nugget, 0, min(vgmdef$gamma))
  init_range <- 0.1*(max(vgmdef$dist)-min(vgmdef$dist))   # 0.10 times the length of the central axis through the area
  init_sill <- mean(c(max(vgmdef$gamma), median(vgmdef$gamma)))
  if(identical(model,''))
    model <- c("Sph","Exp","Gau") # as.character(vgm()$short)

  fitModels <- function(psill, range, nugget, ...) {
    serror <-Inf
    fitmodel <- NULL
    mfit <- NULL
    for(m in model) {
      if(!is.na(fixed.range)) range <- fixed.range
      initSv <- if(!fit.nugget) vgm(psill,m,range,...) else vgm(psill,m,range,nugget,...)

      fit.method <- list(...)$fit.method
      if(is.null(fit.method)) fit.method <- 7
      tryCatch(suppressWarnings(mfit <- fit.variogram(vgmdef, initSv, fit.ranges = is.na(fixed.range), fit.method = fit.method)),
               error=function(e) mfit <<- NULL)
      if(is.null(mfit)) next

      sserr <- attr(mfit, "SSErr")
      if(is.na(sserr)) next

      if(sserr < serror) {
        serror <- attr(mfit, "SSErr")
        fitmodel <- mfit
      }
      # 			if(!attr(mfit, "singular")) break
    }
    # 		return(list(bins=vgmdef, model=fitmodel, sserr=serror))
    return(list(model=fitmodel, sserr=serror))
  }

  # try with initial parameters
  psill <- init_sill-init_nugget
  range <- init_range
  nugget <- init_nugget
  fModel <- fitModels(psill, range, nugget, ...)
  if(!is.null(fModel) && !is.null(fModel$model)) {
    fModel$bins <- vgmdef
    if(fig) print(plot(fModel$bins, fModel$model))
    return(fModel)
  }

  # try with different parameter combinations
  serror <-Inf
  # comb <- c(1/3, 3, 1/9, 9, 1/27, 27, 1/100, 100, 1/500, 500)
  comb <- c(1/3, 3, 1/9, 9, 1/27, 27)
  for(n in nugget*comb) {
    for(r in range*comb) {
      for(p in psill*comb) {
        mfit <- fitModels(p, r, n, ...)
        if(!is.null(attr(mfit, "SSErr")) && attr(mfit, "SSErr") < serror) {
          serror <- attr(mfit, "SSErr")
          fModel <- mfit
          if(!attr(fModel, "singular")) {
            fModel$bins <- vgmdef
            if(fig) plot(fModel$bins, fModel$model)
            return(fModel)
          }
        }
      }
    }
  }

  fModel$bins <- vgmdef
  if(fig) plot(fModel$bins, fModel$model)
  return(fModel)
}


## .calcSvCloudAreaByPointVgm: [internal use only] Variogram scatter for areas according to point scale variogram model. ----
# discretePoints: discretized points of areas, four columns: areaId,x,y,weight.
# ptVgmModel: point-scale semivariogram (gstat vgm).
.calcSvCloudAreaByPointVgm <- function(discretePoints, ptVgmModel, areaDistByCentroid, areaDistByPts) {
  uId <- sort(unique(discretePoints[,1]))

  hasCluster <- !is.null(getOption("ataKrigCluster"))
  if(!hasCluster) {
    dg <- matrix(nrow=length(uId)*(length(uId)-1)/2, ncol=2)
    n <- 1
    for(i in 1:(length(uId)-1)) {
      mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[i]][[i]])
      areaIndexI <- which(discretePoints[,1] == uId[i])
      g11 <- sum(outer(discretePoints[areaIndexI,4], discretePoints[areaIndexI,4]) * mSvar)
      for(j in (i+1):length(uId)) {
        mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[j]][[j]])
        areaIndexJ <- which(discretePoints[,1] == uId[j])
        g22 <- sum(outer(discretePoints[areaIndexJ,4],discretePoints[areaIndexJ,4]) * mSvar)

        mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[i]][[j]])
        g12 <- sum(outer(discretePoints[areaIndexI,4],discretePoints[areaIndexJ,4]) * mSvar)

        g <- g12 - (g11+g22)/2
        dg[n,] <- c(areaDistByCentroid[i,j], g)
        n <- n+1
      }
    }
  } else {
    dg <- foreach(i = 1:(length(uId)-1), .packages = "gstat") %dopar%  {
      mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[i]][[i]])
      areaIndexI <- which(discretePoints[,1] == uId[i])
      g11 <- sum(outer(discretePoints[areaIndexI,4], discretePoints[areaIndexI,4]) * mSvar)

      n <- 1
      dg <- matrix(nrow=length(uId)-i, ncol=2)
      for(j in (i+1):length(uId)) {
        mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[j]][[j]])
        areaIndexJ <- which(discretePoints[,1] == uId[j])
        g22 <- sum(outer(discretePoints[areaIndexJ,4],discretePoints[areaIndexJ,4]) * mSvar)

        mSvar <- variogramLine(ptVgmModel, covariance=FALSE, dist_vector=areaDistByPts[[i]][[j]])
        g12 <- sum(outer(discretePoints[areaIndexI,4],discretePoints[areaIndexJ,4]) * mSvar)

        g <- g12 - (g11+g22)/2
        dg[n,] <- c(areaDistByCentroid[i,j], g)
        n <- n+1
      }
      return(dg)
    }
    dg <- do.call(rbind, dg)
    ataClusterClearObj()
  }

  dg <- as.data.frame(dg)
  names(dg) <- c("dist", "gamma")

  return(dg)
}


## .calcSvCloud: [internal use only] Variogram scatter for points. ----
#      pts: data.frame(ptx,pty,value)
.calcSvCloud <- function(x, y=x, longlat=FALSE) {
  distM <- spDists(as.matrix(x[,1:2]), as.matrix(y[,1:2]), longlat=longlat)
  gammaM <- distM * NA
  for(i in 1:nrow(x)) {
    gammaM[i,] <- (x[i,3]-y[,3])^2 / 2
  }

  if(identical(x, y)) {
    d <- distM[lower.tri(distM)]
    g <- gammaM[lower.tri(gammaM)]
    dg <- data.frame(d=d, g=g)
  } else {
    dg <- data.frame(d=as.vector(distM), g=as.vector(gammaM))
  }
  return(dg)
}


## .calcSvDgGroup: [internal use only] Group variogram scatter (experience variogram). ----
.calcSvDgGroup <- function(dgScatter, ngroup=12, rd=0.66, removeOdd=FALSE) {
	# if(removeOdd) {
	# 	gmn <- mean(dgScatter$g)
	# 	gst <- sd(dgScatter$g)
	# 	indx <- (dgScatter$g > gmn+3*gst) | (dgScatter$g < gmn-3*gst)
	# 	if(sum(indx) > 0) {
	# 		dgScatter$d <- dgScatter$d[-which(indx)]
	# 		dgScatter$g <- dgScatter$g[-which(indx)]
	# 	}
	# }

	if(ngroup > 0) {
# 		dgGroup <- c()
		dmin <- min(dgScatter$d)
		dmax <- max(dgScatter$d)*rd

		dfrom <- dmin
		dinterval <- (dmax-dmin)/ngroup
		dgGroup <- data.frame(np=rep(0,ngroup), dist=rep(0,ngroup),
						gamma=rep(NA,ngroup), dir.hor=rep(0,ngroup),
						dir.ver=rep(0,ngroup), id=rep('var1',ngroup))
		ncnt <- 0
		for(i in 1:ngroup) {
			dto <- dfrom+dinterval
			indx <- (dgScatter$d >= dfrom & dgScatter$d < dto)

			m <- sum(indx)
			if(m >= 1) {
				ncnt <- ncnt+1
				dgGroup[ncnt,] <- data.frame(np=m, dist=mean(dgScatter$d[indx]),
						gamma=mean(dgScatter$g[indx]), dir.hor=0, dir.ver=0, id='var1')
			}
			dfrom <- dto
		}
		if(ncnt < ngroup) dgGroup <- dgGroup[1:ncnt,]
	} else {
		m <- length(dgScatter$d)
		dgGroup <- data.frame(np=rep(1,m), dist=dgScatter$d, gamma=dgScatter$g,
				dir.hor=rep(0,m), dir.ver=rep(0,m), id=rep('var1',m))
	}

	dgGroup$np <- as.numeric(dgGroup$np)

	row.names(dgGroup) <- NULL
	class(dgGroup) <- c("gstatVariogram", "data.frame")
	return(dgGroup)
}


## .calcCrossSvCloudAreaByPointVgm: [internal use only] Cross-variogram scatter for areas according to point scale variogram model. ----
# x, y: discrete points of areas, data.frame(areaId,ptx,pty,weight).
# xPointVgm, yPointVgm, xyPointCrossVgm: point-scale variogram or cross-variogram (gstat vgm).
.calcCrossSvCloudAreaByPointVgm <- function(x, y, xPointVgm, yPointVgm, xyPointCrossVgm,
                                       areaDistByCentroid12, areaDistByPts1, areaDistByPts2, areaDistByPts12) {
  uId1 <- sort(unique(x[,1]))
  uId2 <- sort(unique(y[,1]))

  hasCluster <- !is.null(getOption("ataKrigCluster"))

  if(!hasCluster) {
    dg <- matrix(nrow=length(uId1)*length(uId2), ncol=2)
    n <- 1
    for(i in 1:length(uId1)) {
      mSvar <- variogramLine(xPointVgm, covariance=FALSE, dist_vector=areaDistByPts1[[i]][[i]])
      areaIndexI <- which(x[,1] == uId1[i])
      g11 <- sum(outer(x[areaIndexI,4], x[areaIndexI,4]) * mSvar)
      for(j in 1:length(uId2)) {
        mSvar <- variogramLine(yPointVgm, covariance=FALSE, dist_vector=areaDistByPts2[[j]][[j]])
        areaIndexJ <- which(y[,1] == uId2[j])
        g22 <- sum(outer(y[areaIndexJ,4],y[areaIndexJ,4]) * mSvar)

        mSvar <- variogramLine(xyPointCrossVgm, covariance=FALSE, dist_vector=areaDistByPts12[[i]][[j]])
        g12 <- sum(outer(x[areaIndexI,4],y[areaIndexJ,4]) * mSvar)

        g <- g12 - (g11+g22)/2
        dg[n,] <- c(areaDistByCentroid12[i,j], g)
        n <- n+1
      }
    }
  } else {
    dg <- foreach(i = 1:length(uId1)) %dopar% {
      mSvar <- variogramLine(xPointVgm, covariance=FALSE, dist_vector=areaDistByPts1[[i]][[i]])
      areaIndexI <- which(x[,1] == uId1[i])
      g11 <- sum(outer(x[areaIndexI,4], x[areaIndexI,4]) * mSvar)

      n <- 1
      dg <- matrix(nrow=length(uId2), ncol=2)
      for(j in 1:length(uId2)) {
        mSvar <- variogramLine(yPointVgm, covariance=FALSE, dist_vector=areaDistByPts2[[j]][[j]])
        areaIndexJ <- which(y[,1] == uId2[j])
        g22 <- sum(outer(y[areaIndexJ,4],y[areaIndexJ,4]) * mSvar)

        mSvar <- variogramLine(xyPointCrossVgm, covariance=FALSE, dist_vector=areaDistByPts12[[i]][[j]])
        g12 <- sum(outer(x[areaIndexI,4],y[areaIndexJ,4]) * mSvar)

        g <- g12 - (g11+g22)/2
        dg[n,] <- c(areaDistByCentroid12[i,j], g)
        n <- n+1
      }
      return(dg)
    }
    dg <- do.call(rbind, dg)
    ataClusterClearObj()
  }

  dg <- as.data.frame(dg)
  names(dg) <- c("dist", "gamma")

  return(dg)
}


## .calcAreaCentroid: [internal use only] Area centroid. ----
# discretePoints: discretized points of areas, four columns: areaId,ptx,pty,weight.
.calcAreaCentroid <- function(discretePoints) {
  uId <- sort(unique(discretePoints[,1]))
  xys <- matrix(nrow=length(uId), ncol=2)
  for(i in 1:length(uId)) {
    indx <- (discretePoints[,1] == uId[i])
    ws <- sum(discretePoints[indx,4])
    x <- sum(discretePoints[indx,2] * discretePoints[indx,4]/ws)
    y <- sum(discretePoints[indx,3] * discretePoints[indx,4]/ws)
    xys[i,] <- c(x,y)
  }

  colnames(xys) <- c("centx","centy")
  xys <- as.data.frame(xys)
  xys$areaId <- uId
  xys <- xys[,c("areaId","centx","centy")]
  return(xys)
}


## .crossName: [internal use only] names for cross-variogram model ----
.crossName <- function(id1, id2) {
  if(id1 == id2) {
    return(id1)
  } else {
    id <- sort(c(id1,id2))
    return(paste(id[1], id[2], sep = "."))
  }
}


## .meanDeconvVgms ----
.meanDeconvVgms <- function(dVgms) {
  dVgm <- dVgms[[1]]
  ngroup <- nrow(dVgms[[1]]$experientialAreaVariogram)
  sserr <- dVgms[[1]]$sserr
  for (k in 2:length(dVgms)) {
    dVgm$pointVariogram <- rbind(dVgm$pointVariogram, dVgms[[k]]$pointVariogram)
    dVgm$areaVariogram <- rbind(dVgm$areaVariogram, dVgms[[k]]$areaVariogram)
    dVgm$experientialAreaVariogram <- rbind(dVgm$experientialAreaVariogram, dVgms[[k]]$experientialAreaVariogram)
    dVgm$regularizedAreaVariogram <- rbind(dVgm$regularizedAreaVariogram, dVgms[[k]]$regularizedAreaVariogram)
    ngroup <- c(ngroup, nrow(dVgms[[k]]$experientialAreaVariogram))
    sserr <- c(sserr, dVgms[[k]]$sserr)
  }

  ngroupFreq <- as.data.frame(table(ngroup))
  ngroup <- ngroupFreq$ngroup[which.max(ngroupFreq$Freq)]
  ngroup <- as.integer(levels(ngroup)[ngroup])

  # w <- 1/sserr
  w <- rep(1, length(dVgms))
  w <- w/sum(w)

  if(length(unique(dVgm$pointVariogram$model)) == 1) {
    # m <- colMeans(dVgm$pointVariogram[,-1])
    m <- colSums(dVgm$pointVariogram[,-1] *
                   matrix(rep(w, ncol(dVgm$pointVariogram)),
                          ncol = ncol(dVgm$pointVariogram)))
    dVgm$pointVariogram <- dVgm$pointVariogram[1,]
    dVgm$pointVariogram[,-1] <- m
  } else {
    stop("Model types different!")
  }

  if(length(unique(dVgm$areaVariogram$model)) == 1) {
    # m <- colMeans(dVgm$areaVariogram[,-1])
    m <- colSums(dVgm$areaVariogram[,-1] *
                   matrix(rep(w, ncol(dVgm$areaVariogram)),
                          ncol = ncol(dVgm$areaVariogram)))
    dVgm$areaVariogram <- dVgm$areaVariogram[1,]
    dVgm$areaVariogram[,-1] <- m
  } else {
    stop("Model types different!")
  }

  expAreaVgm <- c()
  drng <- range(dVgm$experientialAreaVariogram$dist)
  dintv <- seq(drng[1], drng[2], length.out = 20)
  dindex <- findInterval(dVgm$experientialAreaVariogram$dist, dintv, rightmost.closed = T)
  for (k in sort(unique(dindex))) {
    cur <- dVgm$experientialAreaVariogram[dindex == k,]
    w <- cur$np/sum(cur$np)
    expAreaVgm <- rbind(expAreaVgm, data.frame(np = sum(cur$np), dist = sum(w*cur$dist), gamma = sum(w*cur$gamma)))
  }
  expAreaVgm$dir.hor <- 0
  expAreaVgm$dir.ver <- 0
  expAreaVgm$id <- "var1"
  class(expAreaVgm) <- class(dVgm$experientialAreaVariogram)
  dVgm$experientialAreaVariogram <- expAreaVgm

  regAreaVgm <- c()
  drng <- range(dVgm$regularizedAreaVariogram$dist)
  dintv <- seq(drng[1], drng[2], length.out = 20)
  dindex <- findInterval(dVgm$regularizedAreaVariogram$dist, dintv, rightmost.closed = T)
  for (k in sort(unique(dindex))) {
    cur <- dVgm$regularizedAreaVariogram[dindex == k,]
    regAreaVgm <- rbind(regAreaVgm, data.frame(dist = mean(cur$dist), gamma = mean(cur$gamma)))
  }
  class(regAreaVgm) <- class(dVgm$regularizedAreaVariogram)
  dVgm$regularizedAreaVariogram <- regAreaVgm

  dVgm$status <- NULL
  dVgm$sserr <- mean(sserr)

  return(dVgm)
}


## extractPointVgm: extract  point-scale variograms from deconvoluted result. ----
extractPointVgm <- function(g) {
  if(!is(g, "ataKrigVgm")) return(NULL)

  if(hasName(g, "pointVariogram")) {
    return(g$pointVariogram)
  } else {
    for(id in names(g)) {
      if(hasName(g[[id]], "pointVariogram"))
        g[[id]] <- g[[id]]$pointVariogram
    }
    return(g)
  }
}
