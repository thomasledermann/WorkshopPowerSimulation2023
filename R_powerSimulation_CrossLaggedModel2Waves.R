#########################################################################################################################
### P o w e r   A n a l y s i s   f o r   a   C r o s s - L a g g e d   M o d e l   a n d   2   T i m e   P o i n t s ###
### Author: Thomas Ledermann                                                                                          ###
### Created: November 16, 2023                                                                                        ###
### Last update: November 16, 2023                                                                                    ###
#########################################################################################################################

sampleSize <- 100
alphaLevel <- .05		# significance level
nsim <- 10.000		# number of iterations
setSeed <- 123

rx1x2 <- .10	# correlation between X1 and X2
rx1y1 <- .30	# correlation between X1 and Y1
rx1y2 <- .20	# correlation between X1 and Y2
rx2y1 <- .20	# correlation between X2 and Y1
rx2y2 <- .30	# correlation between X2 and Y2
ry1y2 <- .30	# correlation between Y1 and Y2

vX1 <- 1		# VAR(X1), can be set to 1
vX2 <- 1		# VAR(X2), can be set to 1
vY1 <- 1		# VAR(Y1), can be set to 1
vY2 <- 1		# VAR(Y2), can be set to 1

# install and load packages
if(!require("lavaan")) install.packages("lavaan")
if(!require("paramtest")) install.packages("paramtest")
if(!require("simsem")) install.packages("simsem")
if(!require("dplyr")) install.packages("dplyr")

library(lavaan)
library(paramtest)
library(simsem)
library(dplyr)

# Covariance matrix
rMat <- matrix(NA, 4, 4)
rMat[lower.tri(rMat, diag = TRUE)] <- c(1, rx1x2, rx1y1, rx1y2, 1, rx2y1, rx2y2, 1, ry1y2, 1)
rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)]
vNames <- c('x1', 'x2', 'y1', 'y2')                    # Variable names
dimnames(rMat) <- list(vNames, vNames)
sds <- sqrt(c(vX1, vX2, vY1, vY2))
covMat <- sds %*% t(sds) * rMat
covMat

## Estimate the parameters using the covariance matrix
CLM <- '	
	y1 ~ a1*x1 + c21*x2
	y2 ~ c12*x1 + a2*x2
	x1 ~~ vx1*x1 + cx*x2
	x2 ~~ vx2*x2
	y1 ~~ vy1*y1 + cy*y2
	y2 ~~ vy2*y2
'

fit <- sem(CLM, sample.cov = covMat, sample.nobs = 100000)
summary(fit, standardized = TRUE, rsquare = TRUE)

ests <- parameterEstimates(fit)
ests

# extract results from the lavaan output
AR1 <- ests[ests$label == 'a1', 'est']
AR2 <- ests[ests$label == 'a2', 'est']
CL12 <- ests[ests$label == 'c12', 'est']
CL21 <- ests[ests$label == 'c21', 'est']
covX <- ests[ests$label == 'cx', 'est']
covY <- ests[ests$label == 'cy', 'est']
varX1 <- ests[ests$label == 'vx1', 'est']
varX2 <- ests[ests$label == 'vx2', 'est']
varY1 <- ests[ests$label == 'vy1', 'est']
varY2 <- ests[ests$label == 'vy2', 'est']

popEst <- cbind.data.frame(AR1, AR2, CL21, CL12, covX, covY, varX1, varX2, varY1, varY2)
popEst

## Power Simulation
models <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$AR1, '*X1
		Y1 ~ ', .$CL21, '*X2
		Y2 ~ ', .$CL12, '*X1
		Y2 ~ ', .$AR2, '*X2
		X1 ~~ ', .$covX, '*X2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ AR1*X1 + CL21*X2
		Y2 ~ CL12*X1 + AR2*X2
		X1 ~~ covX*X2
		Y1 ~~ covY*Y2
	'
	data.frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY,
		gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
models

powerSim <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSim

# Power estimates and other statistics
round(summaryParam(powerSim$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Power Estimates
pEst <- getPower(powerSim$pSimEst[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'))
pEst



## Find Sample Size
sampleSizes <- seq(25, 800, 50)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 10.00	# number of replications

SimNest <- models %>% do({
	SimN <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY,
	SimEstN = list(SimN))
})

Nest <- getPower(SimNest$SimEstN[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'), nVal = lowerN:upperN)
sampleSizeEst <- findPower(Nest, iv = "N", power = 0.80)
sampleSizeEst

## Evaluation
# No parameter names in the fitModel for additional statistics 
modelsNoNames <- popEst %>% rowwise() %>% do({
	genModel <- paste0('
		Y1 ~ ', .$AR1, '*X1
		Y1 ~ ', .$CL21, '*X2
		Y2 ~ ', .$CL12, '*X1
		Y2 ~ ', .$AR2, '*X2
		X1 ~~ ', .$covX, '*X2
		Y1 ~~ ', .$covY, '*Y2
		X1 ~~ ', .$varX1, '*X1
		X2 ~~ ', .$varX2, '*X2
		Y1 ~~ ', .$varY1, '*Y1
		Y2 ~~ ', .$varY2, '*Y2
	')
	fitModel <-'
		Y1 ~ X1 + X2
		Y2 ~ X1 + X2
		X1 ~~ X2
		Y1 ~~ Y2
	'
	data.frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY,
		gen = genModel, fit = fitModel, stringsAsFactors = FALSE)
})
modelsNoNames 

powerSimNoNames <- modelsNoNames %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY,
	pSimEst = list(pSim))
})

# Power Estimates and other statistics
round(summaryParam(powerSimNoNames$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Coverage
CoverN <- getCoverage(powerSimNoNames$pSimEst[[1]], coverParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'))
CoverN

# Power Estimates
powerEst <- getPower(powerSimNoNames$pSimEst[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'))
powerEst



## Skewness: skewness = (0, 0, 3, 3)
distSk <- bindDist(skewness = c(0, 0, 3, 3), kurtosis = c(0, 0, 0, 0))
powerSimSk <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], indDist = distSk, lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSimSk

# Power Estimates
powerEstSk <- getPower(powerSimSk$pSimEst[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'))
powerEstSk

round(summaryParam(powerSimSk$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Find sample size
sampleSizes <- seq(75, 3500, 200)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 10.00	# number of replications

SimNestSk <- models %>% do({
	SimNSk <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], indDist = distSk, lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	SimEstNSk = list(SimNSk))
})

NestSk <- getPower(SimNestSk$SimEstNSk[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'), nVal = lowerN:upperN)
sampleSizeEstSk <- findPower(NestSk, iv = "N", power = 0.80)
sampleSizeEstSk


## Kurtosis: skewness = (0, 0, 0, 0), kurtosis = (0, 0, 21, 21)
distKu <- bindDist(skewness = c(0, 0, 0, 0), kurtosis = c(0, 0, 21, 21))
powerSimKu <- models %>% do({
	pSim <- sim(nRep = nsim, model = .$fit[1], n = sampleSize, generate = .$gen[1], indDist = distKu, lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	pSimEst = list(pSim))
})

# Point estimates of the population model
powerSimKu

# Power Estimates
powerEstKu <- getPower(powerSimKu$pSimEst[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'))
powerEstKu

round(summaryParam(powerSimKu$pSimEst[[1]], detail = TRUE, alpha = 0.05), 3)

# Find sample size
sampleSizes <- seq(25, 800, 100)
lowerN <- min(sampleSizes)
upperN <- max(sampleSizes)
nRep <- 10.00	# number of replications

SimNestKu <- models %>% do({
	SimNKu <- sim(nRep = NULL, model = .$fit[1], n = rep(sampleSizes, nRep), generate = .$gen[1], indDist = distKu, lavaanfun = "sem")
	data_frame(AR1 = .$AR1, AR2 = .$AR2, CL21 = .$CL21, CL12 = .$CL12, covX = .$covX, covY = .$covY, varY1 = .$varY1, varY2 = .$varY2,
	SimEstNKu = list(SimNKu))
})

NestKu <- getPower(SimNestKu$SimEstNKu[[1]], powerParam = c('AR1', 'AR2', 'CL21', 'CL12', 'covX', 'covY'), nVal = lowerN:upperN)
sampleSizeEstKu <- findPower(NestKu, iv = "N", power = 0.80)
sampleSizeEstKu
