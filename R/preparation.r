## atakrig
## Author: Maogui Hu.


## discretize raster to fine resolution points. ----
# Input:
#   x: raster.
#   cellsize: discretized fine grid size.
#   type: "value", "nodata", "all", whether only valid pixels, or only NODATA pixles, or all pixels extracted.
#   psf: PSF type, "equal", "gau", or user defined PSF matrix (normalized).
#   sigma: standard deviation for Gaussian PSF.
# Output: list(areaValues, discretePoints)
#   areaValues: original pixel coordinates and value, data.frame(areaId,centx,centy,value).
#   discretePoints: discretized pixels points and weight, data.frame(areaId,ptx,pty,weight).
discretizeRaster <- function(x, cellsize, type="value", psf="equal", sigma=2) {
  if(!is(x, "RasterLayer")) {
    stop("x isn't raster object!")
  }

  type <- tolower(type)
  if(type == "value") {
    p <- data.frame(rasterToPoints(x))
  } else if(type == "nodata") {
    r0 <- is.na(x)
    p <- data.frame(rasterToPoints(r0, fun = function(x) x > 0))
    p[,3] <- NA
  } else if(type == "all"){
    p1 <- data.frame(rasterToPoints(x))
    names(p1) <- c("x","y","value")

    r0 <- is.na(x)
    p <- data.frame(rasterToPoints(r0, fun = function(x) x > 0))
    p[,3] <- NA
    names(p) <- c("x","y","value")

    p <- rbind(p1, p)
  } else {
    error("invalid type.")
  }

  p <- data.frame(areaId=1:nrow(p), p)

  blkResoX <- xres(x)
  blkResoY <- yres(x)
  xResoNum <- round(blkResoX / cellsize)
  yResoNum <- round(blkResoY / cellsize)

  if(is.numeric(psf)) {
    if(!all(dim(psf) == c(yResoNum, xResoNum))) {
      stop("Incorrect psf dimension!")
    }
    w <- psf
  } else if(psf == "equal") {
    w <- 1/(xResoNum*yResoNum)
  } else if(psf == "gau") {
    xy <- expand.grid(x=(1:xResoNum)-ceiling(xResoNum/2), y=(1:yResoNum)-ceiling(yResoNum/2))
    w <- 1/(2 * pi * sigma^2) * exp(-(xy$x^2 + xy$y^2)/(2 * sigma^2))
    # w <- matrix(w/sum(w), ncol = xResoNum, byrow = TRUE)
    w <- w/sum(w)
  } else {
    stop("Incorrect psf parameter!")
  }

  discretePoints <- lapply(1:nrow(p), function(i) {
    xy <- expand.grid(list(ptx = p$x[i] - 0.5 * blkResoX + (1:xResoNum-0.5) * cellsize,
                           pty = p$y[i] - 0.5 * blkResoY + (1:yResoNum-0.5) * cellsize))
    data.frame(areaId = p$areaId[i], xy, weight=w)
  })
  discretePoints <- do.call(rbind, discretePoints)

  colnames(p)[2:4] <- c("centx","centy","value")
  rslt <- list(areaValues = p, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")

  return(rslt)
}


discretizePolygon <- function(x, cellsize, id=NULL, value=NULL, showProgressBar=FALSE) {

  if(is.null(id)) {
    id <- "_id_"
    x$`_id_` <- 1:length(x)
  }
  if(is.null(value)) {
    value <- "_value_"
    x$`_value_` <- NA
  }

  regularsample <- function(x, cellsize, iter=0) {
    bb <- bbox(x)
    xs <- tryCatch(seq(bb[1,1] + 0.5*cellsize, bb[1,2], by=cellsize), error = function(e) mean(bb[1,]))
    ys <- tryCatch(seq(bb[2,1] + 0.5*cellsize, bb[2,2], by=cellsize), error = function(e) mean(bb[2,]))
    xy <- expand.grid(xs, ys)
    names(xy) <- c("x","y")

    pts <- SpatialPoints(xy, CRS(proj4string(x)))[x]
    if(length(pts) == 0 && iter <= 0) iter <- 10

    if(iter > 0) {
      frac <- 0.5/(iter + 1)
      n <- nrow(xy)
      for (i in 1:iter) {
        offsetxy <- xy + frac * i * data.frame(x=rep(cellsize,n), y=rep(cellsize,n))
        pts0 <- SpatialPoints(offsetxy, CRS(proj4string(x)))[x]
        if(length(pts0) > length(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=rep(cellsize,n), y=-rep(cellsize,n))
        pts0 <- SpatialPoints(offsetxy, CRS(proj4string(x)))[x]
        if(length(pts0) > length(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=-rep(cellsize,n), y=rep(cellsize,n))
        pts0 <- SpatialPoints(offsetxy, CRS(proj4string(x)))[x]
        if(length(pts0) > length(pts)) pts <- pts0

        offsetxy <- xy + frac * i * data.frame(x=-rep(cellsize,n), y=-rep(cellsize,n))
        pts0 <- SpatialPoints(offsetxy, CRS(proj4string(x)))[x]
        if(length(pts0) > length(pts)) pts <- pts0
      }
    }

    return(pts)
  }

  if(class(x@data[,id]) == "factor") {
    x@data[,id] <- as.character(x@data[,id])
  }

  if(showProgressBar) pb <- txtProgressBar(max = length(x), width = 50, style = 3)
  discretePoints <- lapply(1:length(x), function(i) {
    if(showProgressBar) setTxtProgressBar(pb, i)
    pts <- regularsample(x[i,], cellsize = cellsize)
    xys <- cbind(data.frame(areaId=x@data[i,id]), coordinates(pts))
    xys$weight <- 1/nrow(xys)
    return(xys)
  })
  if(showProgressBar) close(pb)
  discretePoints <- do.call(rbind, discretePoints)
  names(discretePoints)[2:3] <- c("ptx", "pty")

  areaValues <- data.frame(areaId=x@data[,id], value=x@data[,value])

  centxy <- rgeos::gCentroid(x, byid = TRUE)
  areaValues <- cbind(areaValues, centxy@coords)
  names(areaValues)[3:4] <- c("centx", "centy")
  areaValues <- areaValues[,c("areaId","centx","centy","value")]

  if(class(areaValues$areaId) == "factor") areaValues$areaId <- as.character(areaValues$areaId)
  if(class(discretePoints$areaId) == "factor") discretePoints$areaId <- as.character(discretePoints$areaId)

  rslt <- list(areaValues = areaValues, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")

  return(rslt)
}


isValidDiscreteAreaObj <- function(x){
  if(!all(sort(names(x))) == c("areaValues", "discretePoints")) return(FALSE)
  if(!all(names(x$areaValues) == c("areaId","centx","centy","value"))) return(FALSE)
  if(!all(names(x$discretePoints) == c("areaId","ptx","pty","weight"))) return(FALSE)
  return(TRUE)
}


updateDiscreteAreaValue <- function(x, newval) {
  for (i in 1:nrow(newval)) {
    indx <- x$areaValues$areaId == newval$areaId
    x$areaValues$value[indx] <- newval$value[i]
  }
  return(x)
}


## select discretized data by area id. ----
subsetDiscreteArea <- function(x, selAreaId, revSel=FALSE) {
  if(!revSel) {
    rslt <- list(areaValues = x$areaValues[x$areaValues[,1] %in% selAreaId,,drop=FALSE],
         discretePoints = x$discretePoints[x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  } else {
    rslt <- list(areaValues = x$areaValues[!x$areaValues[,1] %in% selAreaId,,drop=FALSE],
         discretePoints = x$discretePoints[!x$discretePoints[,1] %in% selAreaId,,drop=FALSE])
  }

  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}


## rbindDiscreteArea ----
# Input:
#   x, y: discretized area, list(areaValues, discretePoints):
#       areaValues: sample values, data.frame(areaId,centx,centy,value).
#       discretePoints: discretized area-samples, data.frame(areaId,ptx,pty,weight), weight is normalized.
rbindDiscreteArea <- function(x, y) {
  if(any(match(unique(y$areaValues$areaId), unique(x$areaValues$areaId)))) {
    x$areaValues$areaId <- paste0("x_", x$areaValues$areaId)
    x$discretePoints$areaId <- paste0("x_", x$discretePoints$areaId)

    y$areaValues$areaId <- paste0("y_", y$areaValues$areaId)
    y$discretePoints$areaId <- paste0("y_", y$discretePoints$areaId)
  }

  xy <- list(areaValues = rbind(x$areaValues, y$areaValues),
             discretePoints = rbind(x$discretePoints, y$discretePoints))
  class(xy) <- c("list", "discreteArea")
  return(xy)
}


## writeValuesBlock: Write values to a file. Modified from raster::writeValues(). ----
# Input:
#   x: RasterLayer.
#   v: matrix. Values to be written.
#   start: Integer. Row and col number (counting starts at 1) from where to start writing v.
# Note: the usage is similar to raster::writeValues(), but only gdal-based driver (such as GTiff) supported!
# Example:
#   r <- raster(system.file("external/test.grd", package="raster"))
#   v <- getValuesBlock(r, row=80, nrows=10, col=10, ncols=20)
#   v <- matrix(v, nrow = 10)
#   s <- raster(r)
#   s <- writeStart(s, filename='test.tif', format='GTiff', overwrite=TRUE)
#   s <- writeValuesBlock(s, v, start=c(80,10))
#   s <- writeStop(s)
writeValuesBlock <- function(x, v, start=c(1,1)) {
  v[is.infinite(v)] <- NA

  datanotation <- x@file@datanotation
  if (substr(datanotation,1,1) != 'F') {
    v <- round(v)
    size <- substr(datanotation,4,4)
    if (substr(datanotation, 1, 3) == 'LOG') {
      v[v != 1] <- 0
    } else if (substr(datanotation, 5, 5) == 'U') {
      v[v < 0] <- NA
      if (size == '1') {
        v[v > 255] <- NA
      } else if (size == '2') {
        v[v > 65535] <- NA
      } else {
        v[v > 4294967295] <- NA
      }
    } else {
      if (size == '1') {
        v[v < -128] <- NA
        v[v > 127] <- NA
      } else if (size == '2') {
        v[v < -32768] <- NA
        v[v > 32767] <- NA
      } else {
        v[v < -2147483648] <- NA
        v[v > 2147483647] <- NA
      }
    }
  }

  rsd <- stats::na.omit(v) # min and max values
  if (length(rsd) > 0) {
    x@data@min <- min(x@data@min, rsd)
    x@data@max <- max(x@data@max, rsd)
  }

  driver <- x@file@driver

  if ( driver == 'gdal' ) {
    v[is.na(v)] <- x@file@nodatavalue
    if(is.matrix(v)) v <- matrix(v, nrow = ncol(v))
    gd <- rgdal::putRasterData(x@file@transient, v, band=1, offset=start-1)
  } else if ( driver == 'big.matrix') {
    b <- attr(x@file, 'big.matrix')
    nrows <- length(v) / ncol(x)
    if(is.vector(v)) {
      nfrom <- nrow(x)*(start[1]-1) + start[2]
      nto <- nfrom+length(v)-1
      nr <- ceiling((nfrom:nto)/nrow(x))
      nc <- (nfrom:nto)%% nrow(x)
      b[cbind(nr, nc)] <- v
    }
    if(is.matrix(v)) {
      b[start[1]+(1:nrow(v))-1, start[2]+(1:ncol(v))-1] <- v
    }
  } else {
    stop("Only gdal-based driver supported!")
  }
  return(x)
}

