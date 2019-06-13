library(atakrig)
library(rgdal)
library(gstat)

rpath <- system.file("extdata", package="rtop")
observations <- readOGR(rpath, "observations")

hist(observations$obs)
observations$obs <- log(observations$obs)

obs.discrete <- discretizePolygon(observations, cellsize=1500, id="ID", value="obs")

pointsv <- deconvPointVgm(obs.discrete, model=c("Sph","Gau","Exp"), ngroup=12, rd=0.75, maxSampleNum = Inf, fig=T)

pred.cv <- ataKriging.cv(obs.discrete, nfold=nrow(obs.discrete$areaValues), pointsv)
pred.cv$pred <- exp(pred.cv$pred)
pred.cv$obs <- exp(pred.cv$value)

# pearson correlation
cor(pred.cv$obs, pred.cv$pred)

# mae
mean(abs(pred.cv$obs - pred.cv$pred))

# rmse
sqrt(mean((pred.cv$obs - pred.cv$pred)^2))


# prediction
predictionLocations <- readOGR(rpath, "predictionLocations")
pred.discrete <- discretizePolygon(predictionLocations, cellsize = 1500, idfield = "ID")

pred <- ataKriging(obs.discrete, pred.discrete, pointsv$pointVariogram, nmax = 10)
pred$pred <- exp(pred$pred)
pred$var <- exp(pred$var)
