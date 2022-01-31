# ----------------------------------------
# Our method
# ----------------------------------------
modFit <- IntCox(svec = seq(5, 30, by=5), rvec = 8:1, dat, mu = .1, quiet = TRUE, silent = FALSE, rho0=100)

valerrsOurs <- array(0, dim=c(length(modFit$s), length(modFit$r),J))
for(kk in 1:(length(modFit$s))){
  for(jj in 1:length(modFit$r)){
    for(ll in 1:J){
      valerrsOurs[kk,jj,ll] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%modFit$beta[,ll,kk,jj])
    }
  }
}

out <- which(apply(valerrsOurs, c(1,2),sum) == min(apply(valerrsOurs, c(1,2),sum)), arr.ind=TRUE)
betaLR <- modFit$beta[,,out[1,1], out[1,2]]


out <- which(apply(valerrsOurs[,2,,drop=F], 1,sum) == min(apply(valerrsOurs[,2,,drop=F], 1,sum)), arr.ind=TRUE)
betaLR2 <- modFit$beta[,,out, 2]
out <- which(apply(valerrsOurs[,4,,drop=F], 1,sum) == min(apply(valerrsOurs[,4,,drop=F], 1,sum)), arr.ind=TRUE)
betaLR4 <- modFit$beta[,,out, 4]
out <- which(apply(valerrsOurs[,6,,drop=F], 1,sum) == min(apply(valerrsOurs[,6,,drop=F], 1,sum)), arr.ind=TRUE)
betaLR6 <- modFit$beta[,,out,6]
out <- which(apply(valerrsOurs[,7,,drop=F], 1,sum) == min(apply(valerrsOurs[,7,,drop=F], 1,sum)), arr.ind=TRUE)
betaLR7 <- modFit$beta[,,out, 7]
out <- which(apply(valerrsOurs[,8,,drop=F], 1,sum) == min(apply(valerrsOurs[,8,,drop=F], 1,sum)), arr.ind=TRUE)
betaLR8 <- modFit$beta[,,out, 8]



# ------------------------------------------------------------------
# Fit separate L2 penalized models 
# ------------------------------------------------------------------
RidgeBeta <- matrix(0, nrow=p, ncol=J)
for(jj in 1:J){
  outObj <- cbind(as.vector(dat[[jj]]$time), as.vector(dat[[jj]]$status))
  colnames(outObj) <- c("time","status")
  temp1 <- glmnet(dat[[jj]]$X, outObj, family = "cox", alpha=0, lambda.min.ratio = 0.001, standardize = FALSE)
  valerrs1 <- rep(0, 100)
  for(kk in 1:length(temp1$lambda)){
    valerrs1[kk] <- coxnet.deviance(y=Surv(datVal[[jj]]$time, datVal[[jj]]$status), pred = datVal[[jj]]$X%*%temp1$beta[,kk])
  }
  cat(jj, "\n")
  RidgeBeta[,jj] <- as.matrix(coefficients(temp1, s=temp1$lambda[which(valerrs1 == min(valerrs1))]))
}

valerrsridge <- matrix(0, J, J)
for(kk in 1:J){
  ridgeBetaLRtemp <- svd2(RidgeBeta)$u[,1:kk, drop=FALSE]%*%diag(svd2(RidgeBeta)$d[1:kk], kk)%*%t(svd2(RidgeBeta)$v[,1:kk, drop=FALSE])
  for(ll in 1:J){
    valerrsridge[kk ,ll ] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%ridgeBetaLRtemp[,ll])
  }
}

# ------------------------------------------------------------------
# Fit separate L1 penalized models
# ------------------------------------------------------------------
LassoBeta <- matrix(0, nrow=p, ncol=J)
for(jj in 1:J){
  outObj <- cbind(as.vector(dat[[jj]]$time), as.vector(dat[[jj]]$status))
  colnames(outObj) <- c("time","status")
  temp1 <- glmnet(dat[[jj]]$X, outObj, family = "cox", alpha=1, lambda.min.ratio = 0.001, standardize = FALSE)
  valerrs1 <- rep(Inf, 100)
  for(kk in 1:length(temp1$lambda)){
    valerrs1[kk] <- coxnet.deviance(y=Surv(datVal[[jj]]$time, datVal[[jj]]$status), pred = datVal[[jj]]$X%*%temp1$beta[,kk])
  }
  cat(jj, "\n")
  LassoBeta[,jj] <- as.matrix(coefficients(temp1, s=temp1$lambda[which(valerrs1 == min(valerrs1))]))
}

valerrslasso <- matrix(0, J, J)
for(kk in 1:J){
  LassoBetaLRtemp <- svd2(LassoBeta)$u[,1:kk, drop=FALSE]%*%diag(svd2(LassoBeta)$d[1:kk], kk)%*%t(svd2(LassoBeta)$v[,1:kk, drop=FALSE])
  for(ll in 1:J){
    valerrslasso[kk ,ll ] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%LassoBetaLRtemp[,ll])
  }
}

# ------------------------------------------------------------------
# Fit separate elastic net penalized models
# ------------------------------------------------------------------
ElasticNetBeta <- matrix(0, nrow=p, ncol=J)
for(jj in 1:J){
  outObj <- cbind(as.vector(dat[[jj]]$time), as.vector(dat[[jj]]$status))
  colnames(outObj) <- c("time","status")
  temp1 <- glmnet(dat[[jj]]$X, outObj, family = "cox", alpha=0.5, lambda.min.ratio = 0.001, standardize = FALSE)
  valerrs1 <- rep(Inf, 100)
  for(kk in 1:length(temp1$lambda)){
    valerrs1[kk] <- coxnet.deviance(y=Surv(datVal[[jj]]$time, datVal[[jj]]$status), pred = datVal[[jj]]$X%*%temp1$beta[,kk])
  }
  cat(jj, "\n")
  ElasticNetBeta[,jj] <- as.matrix(coefficients(temp1, s=temp1$lambda[which(valerrs1 == min(valerrs1))]))
}




# ----------------------------------------------------------------
# Reduced rank fit 
# ----------------------------------------------------------------
modfit_RR <- RRCox.cv(dat, nlambda1 = 100, nlambda2 = 0,
                     delta1 = .1, delta2 = 0, 
                     nfolds = NULL, 
                     standardize.x = FALSE,
                     max.iter = 500,
                     inner.quiet = TRUE, 
                     quiet = FALSE,
                     inner.tol = 1e-7)

valerrsOursRR <- array(0, dim=c(length(modfit_RR$lambda.vec), length(modfit_RR$lambda2.vec),J))
for(kk in 1:length(modfit_RR$lambda.vec)){
  for(jj in 1:length(modfit_RR$lambda2.vec)){
    for(ll in 1:J){
      valerrsOursRR[kk,jj,ll] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%modfit_RR$beta[kk,jj,,ll])
    }
  }
}

out <- which(apply(valerrsOursRR, c(1,2),sum) == min(apply(valerrsOursRR, c(1,2),sum)), arr.ind=TRUE)
betaRR <- modfit_RR$beta[out[1,1], out[1,2],,]

# ----------------------------------------------------------------
# Reduced rank fit 
# ----------------------------------------------------------------
modfit_GL <- RRCox.cv(dat, nlambda1 = 0, nlambda2 = 100,
                     delta1 = 0, delta2 = .1, 
                     nfolds = NULL, 
                     standardize.x = FALSE,
                     max.iter = 500,
                     inner.quiet = TRUE, 
                     quiet = FALSE,
                     inner.tol = 1e-7)

valerrsOursGL <- array(0, dim=c(length(modfit_GL$lambda.vec), length(modfit_GL$lambda2.vec),J))

for(jj in 1:length(modfit_GL$lambda2.vec)){
    for(ll in 1:J){
      valerrsOursGL[1,jj,ll] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%modfit_GL$beta[1,jj,,ll])
  }
}

out <- which(apply(valerrsOursGL, c(1,2),sum) == min(apply(valerrsOursGL, c(1,2),sum)), arr.ind=TRUE)
betaGL <- modfit_GL$beta[out[1,1], out[1,2],,]

# ----------------------------------------------------------------
# Rowsparse fit 
# ----------------------------------------------------------------
modfit_RRGL <- RRCox.cv(dat, nlambda1 = 25, nlambda2 = 25,
                      delta1 = .01, delta2 = 0.1, 
                      nfolds = NULL, 
                      standardize.x = FALSE,
                      max.iter = 1e3,
                      inner.quiet = TRUE, 
                      quiet = FALSE,
                      inner.tol = 1e-7)

valerrsOurs <- array(0, dim=c(length(modfit_RRGL$lambda.vec), length(modfit_RRGL$lambda2.vec),J))
for(kk in 1:(length(modfit_RRGL$lambda.vec))){
  for(jj in 1:length(modfit_RRGL$lambda2.vec)){
    cat(sum(svd(modfit_RRGL$beta[kk,jj,,])$d > 1e-8),"\n")
    for(ll in 1:J){
      valerrsOurs[kk,jj,ll] <- coxnet.deviance(y=Surv(datVal[[ll]]$time, datVal[[ll]]$status), pred = datVal[[ll]]$X%*%modfit_RRGL$beta[kk,jj,,ll])
    }
  }
}

out <- which(apply(valerrsOurs, c(1,2),sum) == min(apply(valerrsOurs, c(1,2),sum)), arr.ind=TRUE)
betaRRGL <- modfit_RRGL$beta[out[1,1], out[1,2],,]


store <- list(
  "RidgeBeta" = RidgeBeta,
  "LassoBeta" = LassoBeta,
  "ElasticNetBeta" = ElasticNetBeta,
  "betaRR" = betaRR,
  "betaGL" = betaGL,
  "betaRRGL" = betaRRGL,
  "valerrsOurs" = valerrsOurs,
  "valerrsOursGL" = valerrsOursGL,
  "valerrsOursRR" = valerrsOursRR,
  "betaLR" = betaLR,
  "betaLR2" = betaLR2,
  "betaLR4" = betaLR4,
  "betaLR6" = betaLR6,
  "betaLR7" = betaLR7,
  "betaLR8" = betaLR8,
  "SigmaX" = SigmaX,
  "valerrsOurs" = valerrsOurs,
  "beta" = beta
)

saveRDS(store, paste(savefile, uu,".RDS", sep=""))


# ---------------------------------------------------------------------------
# Compute statistics and save as separate file
# ---------------------------------------------------------------------------
jours_Err <- rep(0, 6 + J)
sours_Err <- rep(0, 6 + J)
lrours_Err <- rep(0, 6 + J)
l2_Err <- rep(0, 6 + J)
l1_Err <- rep(0, 6 + J)
l2proj_Err <- rep(0, 6 + J)
l1proj_Err <- rep(0, 6 + J)
NCours_Err <- rep(0, 6 + J)
NCours2_Err <- rep(0,6 + J)
NCours4_Err <- rep(0,6 + J)
NCours6_Err <- rep(0,6 + J)
NCours7_Err <- rep(0,6 + J)
NCours8_Err <- rep(0,6 + J)




library(MASS)
library(survival)
library(hdnom)
ridgeBeta <- as.matrix(store$RidgeBeta)
ridgeBetaLR <- svd2(ridgeBeta)$u[,1:r, drop=FALSE]%*%diag(svd2(ridgeBeta)$d[1:r], r)%*%t(svd2(ridgeBeta)$v[,1:r, drop=FALSE])
lassoBeta <- as.matrix(store$LassoBeta)
lassoBetaLR <- svd2(lassoBeta)$u[,1:r, drop=FALSE]%*%diag(svd2(lassoBeta)$d[1:r], r)%*%t(svd2(lassoBeta)$v[,1:r, drop=FALSE])
lrBeta <- store$betaRR
sBeta <- store$betaGL
jBeta <- store$betaRRGL
beta <- store$beta
ncBeta <- store$betaLR
ncBeta2 <- store$betaLR2
ncBeta4 <- store$betaLR4
ncBeta6 <- store$betaLR6
ncBeta7 <- store$betaLR7
ncBeta8 <- store$betaLR8


ME <- function(betaInput){
  sum(diag(tcrossprod(crossprod(beta - betaInput, SigmaX), t(beta - betaInput))))
}

l2_Err[1] <- ME(ridgeBeta)
l2proj_Err[1] <- ME(ridgeBetaLR)
l1_Err[1] <- ME(lassoBeta)
l1proj_Err[1] <- ME(lassoBetaLR)
lrours_Err[1] <- ME(lrBeta)
sours_Err[1] <- ME(sBeta)
jours_Err[1] <- ME(jBeta)
NCours_Err[1] <- ME(ncBeta)
NCours2_Err[1] <- ME(ncBeta2)
NCours4_Err[1] <- ME(ncBeta4)
NCours6_Err[1] <- ME(ncBeta6)
NCours7_Err[1] <- ME(ncBeta7)
NCours8_Err[1] <- ME(ncBeta8)


Conc <- function(betaInput){
  mean(sapply(1:3, function(x){summary(coxph(Surv(datTest[[x]]$Tlat, rep(1, length(datTest[[x]]$Tlat))) ~ datTest[[x]]$X%*%betaInput[,x]))[14]$concordance[1]}))
}

l2_Err[2] <- Conc(ridgeBeta)
l2proj_Err[2] <- Conc(ridgeBetaLR)
l1_Err[2] <- Conc(lassoBeta)
l1proj_Err[2] <- Conc(lassoBetaLR)
lrours_Err[2] <- Conc(lrBeta)
sours_Err[2] <- Conc(sBeta)
jours_Err[2] <- Conc(jBeta)
NCours_Err[2] <- Conc(ncBeta)
NCours2_Err[2] <- Conc(ncBeta2)
NCours4_Err[2] <- Conc(ncBeta4)
NCours6_Err[2] <- Conc(ncBeta6)
NCours7_Err[2] <- Conc(ncBeta7)
NCours8_Err[2] <- Conc(ncBeta8)


getBrierMedian <- function(betaIn){

  outs <- list(NA)
  for(jj in 1:J){
    MedTime <- median(datTest[[jj]]$Tlat)
    temp2 <-  (glmnet_basesurv(dat[[jj]]$time, dat[[jj]]$status, lp = dat[[jj]]$X%*%betaIn[,jj], times.eval = MedTime, centered = FALSE)$base_surv)
    expXb <- exp(datTest[[jj]]$X%*%betaIn[,jj])
    survProb <- Reduce(cbind, lapply(temp2, function(x){x^expXb}))
    tempTime <- datTest[[jj]]$Tlat
    keep1 <- (1*(tempTime > MedTime) - survProb)^2
    outs[[jj]] <- keep1
  }
  return(sum(unlist(lapply(outs, mean)))/J)

}

getBrier75 <- function(betaIn){

  outs <- list(NA)
  for(jj in 1:J){
    qTime <- quantile(datTest[[jj]]$Tlat, .75)
    temp2 <-  (glmnet_basesurv(dat[[jj]]$time, dat[[jj]]$status, lp = dat[[jj]]$X%*%betaIn[,jj], times.eval = qTime, centered = FALSE)$base_surv)
    expXb <- exp(datTest[[jj]]$X%*%betaIn[,jj])
    survProb <- Reduce(cbind, lapply(temp2, function(x){x^expXb}))
    tempTime <- datTest[[jj]]$Tlat
    keep1 <- (1*(tempTime > qTime) - survProb)^2
    outs[[jj]] <- keep1
  }
  return(sum(unlist(lapply(outs, mean)))/J)
  
}

getBrier25 <- function(betaIn){

  outs <- list(NA)
  for(jj in 1:J){
    qTime <- quantile(datTest[[jj]]$Tlat, .25)
    temp2 <-  (glmnet_basesurv(dat[[jj]]$time, dat[[jj]]$status, lp = dat[[jj]]$X%*%betaIn[,jj], times.eval = qTime, centered = FALSE)$base_surv)
    expXb <- exp(datTest[[jj]]$X%*%betaIn[,jj])
    survProb <- Reduce(cbind, lapply(temp2, function(x){x^expXb}))
    tempTime <- datTest[[jj]]$Tlat
    keep1 <- (1*(tempTime > qTime) - survProb)^2
    outs[[jj]] <- keep1
  }
  return(sum(unlist(lapply(outs, mean)))/J)
  
}

l2_Err[3] <- getBrierMedian(ridgeBeta)
l2proj_Err[3] <- getBrierMedian(ridgeBetaLR)
l1_Err[3] <- getBrierMedian(lassoBeta)
l1proj_Err[3] <- getBrierMedian(lassoBetaLR)
lrours_Err[3] <- getBrierMedian(lrBeta)
sours_Err[3] <- getBrierMedian(sBeta)
jours_Err[3] <- getBrierMedian(jBeta)
NCours_Err[3] <- getBrierMedian(ncBeta)
NCours2_Err[3] <- getBrierMedian(ncBeta2)
NCours4_Err[3] <- getBrierMedian(ncBeta4)
NCours6_Err[3] <- getBrierMedian(ncBeta6)
NCours7_Err[3] <- getBrierMedian(ncBeta7)
NCours8_Err[3] <- getBrierMedian(ncBeta8)



l2_Err[4] <- getBrier75(ridgeBeta)
l2proj_Err[4] <- getBrier75(ridgeBetaLR)
l1_Err[4] <- getBrier75(lassoBeta)
l1proj_Err[4] <- getBrier75(lassoBetaLR)
lrours_Err[4] <- getBrier75(lrBeta)
sours_Err[4] <- getBrier75(sBeta)
jours_Err[4] <- getBrier75(jBeta)
NCours_Err[4] <- getBrier75(ncBeta)
NCours2_Err[4] <- getBrier75(ncBeta2)
NCours4_Err[4] <- getBrier75(ncBeta4)
NCours6_Err[4] <- getBrier75(ncBeta6)
NCours7_Err[4] <- getBrier75(ncBeta7)
NCours8_Err[4] <- getBrier75(ncBeta8)


l2_Err[5] <- getBrier25(ridgeBeta)
l2proj_Err[5] <- getBrier25(ridgeBetaLR)
l1_Err[5] <- getBrier25(lassoBeta)
l1proj_Err[5] <- getBrier25(lassoBetaLR)
lrours_Err[5] <- getBrier25(lrBeta)
sours_Err[5] <- getBrier25(sBeta)
jours_Err[5] <- getBrier25(jBeta)
NCours_Err[5] <- getBrier25(ncBeta)
NCours2_Err[5] <- getBrier25(ncBeta2)
NCours4_Err[5] <- getBrier25(ncBeta4)
NCours6_Err[5] <- getBrier25(ncBeta6)
NCours7_Err[5] <- getBrier25(ncBeta7)
NCours8_Err[5] <- getBrier25(ncBeta8)



MSE <- function(betaInput){
  mean((beta - betaInput)^2)
}

l2_Err[6] <- MSE(ridgeBeta)
l2proj_Err[6] <- MSE(ridgeBetaLR)
l1_Err[6] <- MSE(lassoBeta)
l1proj_Err[6] <- MSE(lassoBetaLR)
lrours_Err[6] <- MSE(lrBeta)
sours_Err[6] <- MSE(sBeta)
jours_Err[6] <- MSE(jBeta)
NCours_Err[6] <- MSE(ncBeta)
NCours2_Err[6] <- MSE(ncBeta2)
NCours4_Err[6] <- MSE(ncBeta4)
NCours6_Err[6] <- MSE(ncBeta6)
NCours7_Err[6] <- MSE(ncBeta7)
NCours8_Err[6] <- MSE(ncBeta8)


for(j in 1:J){
  MSE <- function(betaInput){crossprod(beta[,j]-  betaInput[,j], SigmaX)%*%(beta[,j]-  betaInput[,j])}
  l2_Err[6+j] <- MSE(ridgeBeta)
  l2proj_Err[6+j] <- MSE(ridgeBetaLR)
  l1_Err[6+j] <- MSE(lassoBeta)
  l1proj_Err[6+j] <- MSE(lassoBetaLR)
  lrours_Err[6+j] <- MSE(lrBeta)
  sours_Err[6+j] <- MSE(sBeta)
  jours_Err[6+j] <- MSE(jBeta)
  NCours_Err[6+j] <- MSE(ncBeta)
  NCours2_Err[6+j] <- MSE(ncBeta2)
  NCours4_Err[6+j] <- MSE(ncBeta4)
  NCours6_Err[6+j] <- MSE(ncBeta6)
  NCours7_Err[6+j] <- MSE(ncBeta7)
  NCours8_Err[6+j] <- MSE(ncBeta8)
}

results <- list(
  "jours_Err" = jours_Err,
  "sours_Err" = sours_Err,
  "lrours_Err" = lrours_Err,
  "l2_Err" = l2_Err,
  "l1_Err" = l1_Err,
  "l2proj_Err" = l2proj_Err,
  "l1proj_Err" = l1proj_Err,
  "NCours_Err" = NCours_Err,
  "NCours2_Err" = NCours2_Err,
  "NCours4_Err" = NCours4_Err,
  "NCours6_Err" = NCours6_Err,
  "NCours7_Err" = NCours7_Err,
  "NCours8_Err" = NCours8_Err
  )

saveRDS(results, file=paste(savefile, uu,"_Summary.RDS", sep=""))


