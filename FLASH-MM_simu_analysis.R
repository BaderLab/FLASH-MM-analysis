##################################################

library(MASS)
library(Matrix)
library(lme4)
library(lmerTest)
library(nebula)
library(FLASHMM)
	



##################################################
##Real data
##Quality control (Cell and gene filtering)

##PBMC data
datafile <- "PBCMdata_ExperimentHubKang_muscat.RData"
load(file = datafile)

##number of celss
nCells <- rowSums(counts > 0)

##number of cells in a group_id
nCellsgrp <- do.call(cbind, 
		tapply(1:ncol(counts), as.factor(coldata$group_id), 
		function(j) rowSums(counts[, j, drop = F] > 0))
		)

dim(nCellsgrp)
head(nCellsgrp)

##number of cells in a cluster_id
nCellscls <- do.call(cbind, 
		tapply(1:ncol(counts), as.factor(coldata$cluster_id), 
		function(j) rowSums(counts[, j, drop = F] > 0))
		)

dim(nCellscls)
		
##number of counts
nCounts <- rowSums(counts)

maxCounts <- 2^20
sum(nCounts > maxCounts)

##nebula filtering
##Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.
cpc <- rowSums(counts)/ncol(counts)
sum(cpc <= 0.005) 

##Filtering
minCells <- 32 
minCellsgrp <- 10
minCellscls <- 5
maxCounts <- 2^20
mincpc <- 0.005

index <- (nCells >= minCells) & (rowSums(nCellsgrp >= minCellsgrp) >= ncol(nCellsgrp)) 
index <- index & (rowSums(nCellscls >= minCellscls) >= 2) #ncol(nCellscls))
index <- index & (nCounts <= maxCounts) 
index <- index & (cpc > mincpc)
sum(index)

counts <- counts[index, ]

rm(index, nCells, nCellsgrp, nCellscls, nCounts, cpc)


##################################################
##Differential analysis
Sys.time()


##parameters
ntrt <- 2  ##number of treatments
ncls <- 12 ##number of clusters (cell-types)
nsam <- 25 ##number of samples (subjects or individuals)
nDEgenes <- 40*ncls  ##number of DE genes specific to a cluster
Ng <- 6000

minbeta <- 0.25 #0.5
maxbeta <- 1
#direct <- "up"
#minbeta <- -1
#maxbeta <- -0.25
#direct <- "dn"
direct <- ""
one.direction <- FALSE #TRUE

var.random <- 0.16

#SEED <- 240305
SEED <- 240403

NBMM.par <- list(ntrt = ntrt, ncls = ncls, nsam = nsam, 
	nGenes = Ng, nDEgenes = nDEgenes, minbeta = minbeta, maxbeta = maxbeta,
	one.direction = one.direction, var.randomeffect = var.random, seed = SEED)

NcList <- c(2e4, 4e4, 6e4, 8e4, 10e4, 12e4)


##########
##lmerTest::lmer with default setting
##lmer p-values
##lmer with default setting


for (Nc in NcList){

NBMM.par$nCells <- Nc
		
#####
##Simulate data
set.seed(SEED)

t1 <- Sys.time()
simudata <- simuRNAseq(counts, nGenes = Ng, nCells = Nc, nsam = nsam, ncls = ncls, ntrt = 2, nDEgenes = nDEgenes, minbeta = minbeta, maxbeta = maxbeta, one.direction = one.direction, var.randomeffect = var.random)
t2 <- Sys.time()

##log transformed counts
Y <- log(1 + simudata$counts)
dim(Y)

##
dat <- simudata$metadata
dat$totalCounts <- colSums(simudata$counts)

dim(dat)
head(dat)
#hist(log(dat$totalCounts))

##
simudata$mean.dispersion <- NULL
simudata$counts <- NULL

#####
flme <- NULL
tlme <- NULL
#plme <- NULL
slme <- NULL
rtlme <- 0

convlme <- NULL
singlme <- NULL

model.formula <- formula(y ~ log(totalCounts) + cls + cls:trt + (1 | sam))
NBMM.par$model.formula <- model.formula

timeUnits <- "secs"
NBMM.par$timeUnits <- timeUnits

for (j in 1:nrow(Y)) {

##the observations
dat$y <- Y[j, ]

##lmer with default setting
##LMM fitting: lme4::lmer 
t1 <- Sys.time()
#lmerFit <- lme4::lmer(y ~ log(totalCounts) + cls + cls:trt + (1 | sam), data = dat)
lmerFit <- lme4::lmer(model.formula, data = dat)
	t2 <- Sys.time()
	rtlme <- rtlme + difftime(t2, t1, units = timeUnits) 
	rm(t1, t2)

#lmerFit@optinfo
#convlme <- c(convlme, length(lmerFit@optinfo$conv$lme4))
convlme <- c(convlme, lmerFit@optinfo$conv$opt)
singlme <- c(singlme, isSingular(lmerFit))

sfit <- summary(lmerFit)$coefficients
flme <- cbind(flme, sfit[, "Estimate"])
tlme <- cbind(tlme, sfit[, "t value"])
slme <- cbind(slme, as.data.frame(VarCorr(lmerFit))$vcov)

rm(sfit, lmerFit)

}##Y

NBMM.par$metadata <- simudata$metadata
NBMM.par$DEgenes <- simudata$DEgenes
NBMM.par$treatment <- simudata$treatment
save(NBMM.par, flme, tlme, slme, rtlme, convlme, singlme,
	file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, "_lmer.RData"))
	
rm(simudata, Y, dat)
	
}##Nc


##########
##lmm

epsilon <- 1e-8
NBMM.par$lmm.epsilon <- epsilon
timeUnits <- "secs"
NBMM.par$timeUnits <- timeUnits

for (Nc in NcList){

NBMM.par$nCells <- Nc
		
#####
##Simulate data
set.seed(SEED)

t1 <- Sys.time()
simudata <- simuRNAseq(counts, nGenes = Ng, nCells = Nc, nsam = nsam, ncls = ncls, ntrt = 2, nDEgenes = nDEgenes, minbeta = minbeta, maxbeta = maxbeta, one.direction = one.direction, var.randomeffect = var.random)
t2 <- Sys.time()
#difftime(t2, t1)

##log transformed counts
Y <- log(1 + simudata$counts)

##total counts per cell (library sizes)
totalCounts <- colSums(simudata$counts)

##
simudata$mean.dispersion <- NULL
simudata$counts <- NULL

##design matrix for fixed effects
X <- model.matrix(~ log(totalCounts) + cls + cls:trt, data = simudata$metadata)

##design matrix for random effects
Z <- model.matrix(~ 0 + sam, data = simudata$metadata)
d <- ncol(Z)

##LMM fit
t1 <- Sys.time()
fit <- lmmfit(Y, X, Z, d=d, epsilon = epsilon, output.cov = F)
	t2 <- Sys.time()
	difftime(t2, t1)
	
rtlmm <- difftime(t2, t1, units = timeUnits)


#####
n <- nrow(X)
XY <- t(Y%*%X)
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y)
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

#XXinv <- ginv(t(X)%*%X)
XXinv <- chol2inv(chol(t(X)%*%X))
Ynorm <- rowSums(Y*Y) #colSums(Y*Y)

rm(X, Y, Z)

t1 <- Sys.time()
fitss <- lmm(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, Ynorm = Ynorm, n = n, d = d, epsilon = epsilon, output.cov = F)
	t2 <- Sys.time()
	
rtlmmSS <- difftime(t2, t1, units = timeUnits)

print(identical(fit, fitss))

NBMM.par$metadata <- simudata$metadata
NBMM.par$DEgenes <- simudata$DEgenes
NBMM.par$treatment <- simudata$treatment
save(NBMM.par, fit, rtlmm, fitss, rtlmmSS,
	file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, "_lmm.RData"))
	
rm(simudata, fit, fitss)
}##Nc

##########
##nebula
#citation("nebula")
#NcList <- c(4e4, 6e4, 8e4, 1e5)

timeUnits <- "secs"
NBMM.par$timeUnits <- timeUnits

for (Nc in NcList){

NBMM.par$nCells <- Nc
		
#####
##Simulate data
set.seed(SEED)

t1 <- Sys.time()
simudata <- simuRNAseq(counts, nGenes = Ng, nCells = Nc, nsam = nsam, ncls = ncls, ntrt = 2, nDEgenes = nDEgenes, minbeta = minbeta, maxbeta = maxbeta, one.direction = one.direction, var.randomeffect = var.random)
t2 <- Sys.time()

##
Y <- simudata$counts

simudata$mean.dispersion <- NULL
simudata$counts <- NULL

##design matrix for fixed effects
totalCounts <- colSums(Y)
X <- model.matrix(~ log(totalCounts) + cls + cls:trt, data = simudata$metadata)

##random effect (sample groups)
Z <- simudata$metadata$sam
table(Z)
length(Z) 

##
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model 
##(the same model as that in the lme4 package).

#model <- "NBGMM"
model <- "NBLMM"
method <- "LN"

NBMM.par$nebula.model <- model
NBMM.par$nebula.method <- method

t1 <- Sys.time()
negbn <- NULL
#The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = Y[, o], id = as.factor(Z[o]), pred = X[o, ], 
	model = model, method = method, cpc = 0.001, covariance = F, output_re = F)
	t2 <- Sys.time()

#difftime(t2, t1)
#timeUnits <- "secs"
rtnebula <- difftime(t2, t1, units = timeUnits)

##
#table(negbn$convergence)

##
NBMM.par$metadata <- simudata$metadata
NBMM.par$DEgenes <- simudata$DEgenes
NBMM.par$treatment <- simudata$treatment

save(NBMM.par, negbn, rtnebula,
	file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, "_nebula.RData"))
	
rm(simudata, Y, X, Z, negbn)
}##Nc


table(negbn$convergence)
# -40  -10    1 
#   1  148 5851 
##
# 1:   The convergence is reached due to a sufficiently small improvement of the function value.
#-10: The convergence is reached because the gradients are close to zero (i.e., the critical point) and no improvement of the function value can be found.
# -30: The convergence fails because the likelihood function returns NaN.
# -40: The convergence fails because the critical point is not reached and no improvement of the function value can be found.

##################################################


	