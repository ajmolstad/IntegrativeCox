# --------------------------------------
# Simulations for RR-cox
# --------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source("~/blue/LRCox/Simulations_R1/Functions/RRCox_PPG.R")
source("~/blue/LRCox/Simulations_R1/Functions/LRCox.R")
soureCpp("~/blue/LRCox/Simulations_R1/Functions/updateBeta.cpp")

library(Matrix)
library(glmnet)
library(MASS)
library(survival)
set.seed(uu)

J <- 12
Ns <- rep(c(1250, 1350, 1450), 4)
nreps <- 100

temp <- expand.grid(p = c(rep(seq(100, 500, by=100), each = nreps)), 
  r = 3)
p <- temp[uu,1]
r <- temp[uu,2]
temp <- NULL
quan.cens <- .35
kappa <- seq(2000, 2110, by = 10)


genCoxDat <- function(uu, p, Ns, quan.cens){
  
  set.seed(uu)
  X <- list(NA)
  SigmaX <- matrix(0, p, p)
  for(j in 1:p){
    for(k in 1:p){
      SigmaX[j,k] <- .7^(abs(j-k))
    }
  }
  
  for(j in 1:J){
    X[[j]] <- mvrnorm(n = Ns[j], mu = rep(0, p), 
                      SigmaX, tol = 1e-06, empirical = FALSE)
  }
  
  cat(dim(X[[1]]), "\n")
  
  simulGomp <- function(datIndex, lambda, alpha, beta, rateC){
    
    N <- dim(X[[datIndex]])[1]
    x <- as.matrix(X[[datIndex]])
    v <- runif(n=N)
    Tlat <- (1/alpha)*log(1 - (alpha*log(v))/(lambda*exp(x%*%beta)))
    
    # censoring times
    if(datIndex%%3 == 0){  
      temp <- quantile(Tlat, quan.cens+.2)
    } else {
      temp <- quantile(Tlat, quan.cens)
    }  
    
    C <- rexp(n=N, rate=1/temp)
    
    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- 1*(Tlat <= C)
    
    # data set
    list("id"=1:N, "time"=time, "status"=status, "X"=x, "Tlat" = Tlat, "C" = C, "linPred" = x%*%beta)
  }
  
  
  get_beta <- function(sigma.temp){
    beta <- matrix(0, nrow=p, ncol=J)
    nonzeroes <- sample(1:p, 20, replace=FALSE)
    temp <- svd(matrix(rnorm(r*J), nrow=r, ncol=J))$v*(sqrt(2)/sqrt(r))
    beta[nonzeroes, ] <- matrix(runif(20*r, 1, 2)*
                                  sample(c(-1, 1), 20*r, replace=TRUE), nrow=20)%*%t(temp)
    return(beta)
  }
  sigma.temp <- NULL
  beta <- get_beta(sigma.temp)
  
  dat <- list(NA)
  for(kk in 1:J){
    alpha <- pi/(600*sqrt(6))
    lambda <- alpha*exp(- 0.5772 - alpha*kappa[kk])
    dat[[kk]] <- simulGomp(datIndex = kk, lambda=lambda, alpha=alpha, beta=beta[,kk], rateC=0.1)
  }
  
  return(list("dat" = dat, "beta" = beta, "SigmaX" = SigmaX))
  
}

simDat <- genCoxDat(uu = uu, p = p, Ns = Ns, quan.cens = quan.cens)
beta <- simDat$beta
dat <- simDat$dat
SigmaX <- simDat$SigmaX
simDat <- NULL

# ----------------------------------------------------
# Split into training, validation, and testing sets
# ----------------------------------------------------
datVal <- list(NA)
datTest <- list(NA)

for(j in 1:J){
  
  ValInds <- (100*((j-1)%%3 + 1) + 1):(100*((j-1)%%3 + 1) + 150)
  TestInds <- (100*((j-1)%%3 + 1) + 151):Ns[j]
  datVal[[j]] <- list(
    "X" = dat[[j]]$X[ValInds,],
    "time" = dat[[j]]$time[ValInds],
    "status" = dat[[j]]$status[ValInds],
    "Tlat" = dat[[j]]$Tlat[ValInds],
    "linPred" = dat[[j]]$linPred[ValInds])
  
  datTest[[j]] <- list(
    "X" = dat[[j]]$X[TestInds,],
    "time" = dat[[j]]$time[TestInds],
    "status" = dat[[j]]$status[TestInds],
    "Tlat" = dat[[j]]$Tlat[TestInds],
    "linPred" = dat[[j]]$linPred[TestInds])
  
  dat[[j]]$X <- dat[[j]]$X[-c(TestInds, ValInds),]
  dat[[j]]$time <- dat[[j]]$time[-c(TestInds, ValInds)]
  dat[[j]]$status <- dat[[j]]$status[-c(TestInds, ValInds)]
  
}

# --------------
# Check R^2
# --------------
# 1 - exp(-(coxnet.deviance(y=Surv(dat[[1]]$time, dat[[1]]$status), pred = rep(0, dim(dat[[1]]$linPred)[1])) - coxnet.deviance(y=Surv(dat[[1]]$time, dat[[1]]$status), pred = dat[[1]]$linPred))/100)
# 1 - exp(-(coxnet.deviance(y=Surv(dat[[2]]$time, dat[[2]]$status), pred = rep(0, dim(dat[[2]]$linPred)[1])) - coxnet.deviance(y=Surv(dat[[2]]$time, dat[[2]]$status), pred = dat[[2]]$linPred))/200)
# 1 - exp(-(coxnet.deviance(y=Surv(dat[[3]]$time, dat[[3]]$status), pred = rep(0, dim(dat[[3]]$linPred)[1])) - coxnet.deviance(y=Surv(dat[[3]]$time, dat[[3]]$status), pred = dat[[3]]$linPred))/300)

savefile <- paste("~/blue/LRCox/Simulations_R1/Model3/Results/Model3_p",p,"_r",r,"_", sep="")
source("~/blue/LRCox/Simulations_R1/Fit_Main.R")

