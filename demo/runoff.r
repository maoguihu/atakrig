library(atakrig)
library(rgdal)

## load demo data from rtop package ----
if (!require("rtop", quietly = TRUE)) message("rtop library is required for demo data.")
rpath <- system.file("extdata", package="rtop")
observations <- readOGR(rpath, "observations")

observations$obs <- observations$QSUMMER_OB/observations$AREASQKM

## point-scale variogram ----
obs.discrete <- discretizePolygon(observations, cellsize=1500, id="ID", value="obs")
pointsv <- deconvPointVgm(obs.discrete, model="Exp", ngroup=12, rd=0.75, fig=TRUE)

## cross validation ----
pred.cv <- ataKriging.cv(obs.discrete, nfold=length(observations), pointsv, showProgress = TRUE)
names(pred.cv)[6] <- "obs"

summary(pred.cv[,c("obs","pred","var")])
cor(pred.cv$obs, pred.cv$pred)            # Pearson correlation
mean(abs(pred.cv$obs - pred.cv$pred))     # MAE
sqrt(mean((pred.cv$obs - pred.cv$pred)^2))# RMSE

## prediction ----
predictionLocations <- readOGR(rpath, "predictionLocations")
pred.discrete <- discretizePolygon(predictionLocations, cellsize = 1500, id = "ID")
pred <- ataKriging(obs.discrete, pred.discrete, pointsv$pointVariogram)
