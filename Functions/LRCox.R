# --------------------------------------------
# Functions for penalty method LRCox
# --------------------------------------------
library(Matrix)
library(glmnet)
library(MASS)
library(survival)
library(Rfast)
library(Rcpp)
library(RcppArmadillo)

IntCox <- function(svec, rvec, dat, mu, quiet, silent, rho0) {
    
  # ------------------------------------------
  # Preliminaries
  # ------------------------------------------
  J <- length(dat)
  q <- dim(dat[[1]]$z)[2]
  
  taus <- list(NA)
  Es <- list(NA)
  Rs <- list(NA)
  Ds <- list(NA)
  SDs <- list(NA)
  
  for (jj in 1:J) {
    sortedInds <- sort(dat[[jj]]$time, decreasing=FALSE, index=TRUE)$ix
    dat[[jj]]$time <- dat[[jj]]$time[sortedInds]
    dat[[jj]]$status <- dat[[jj]]$status[sortedInds]
    dat[[jj]]$X <- dat[[jj]]$X[sortedInds,]
  }
  
  # -----------------------------------------------------
  # Discard censored subjects before first failure time
  # -----------------------------------------------------
  for (jj in 1:J) {
    rmInds <- which(dat[[jj]]$time < min(dat[[jj]]$time[dat[[jj]]$status==1]))
    if (length(rmInds) > 0) {
      dat[[jj]]$time <- dat[[jj]]$time[-rmInds]
      dat[[jj]]$status <- dat[[jj]]$status[-rmInds]
      dat[[jj]]$X <- dat[[jj]]$X[-rmInds,]
    }
  }
  
  
  for (jj in 1:J) {
    taus[[jj]] <- dat[[jj]]$time[which(dat[[jj]]$status==1)]
    Es[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time == x)})
    Ds[[jj]] <- sapply(taus[[jj]], function(x){length(which(dat[[jj]]$time == x))})
    Rs[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time >= x)})
  }
  
  for (jj in 1:J) {
    SDs[[jj]] <- list(NA)
    for (kk in 2:length(taus[[jj]])) {
      SDs[[jj]][[kk-1]] <- setdiff(Rs[[jj]][[kk-1]], Rs[[jj]][[kk]])
    }
  }
  
  
  # ------------------------------------------------------
  # Remove genes with no variability in at least one fold
  # ------------------------------------------------------
  p <- dim(dat[[1]]$X)[2]
  for (jj in 1:J) {
    dat[[jj]]$x <- dat[[jj]]$X
  }
  
  
  getGrads <-  function(dat, beta, taus, Ds, Es, Rs, SDs) {
    J <- length(dat)
    ns <- unlist(lapply(dat, function(x){dim(x$x)[1]}))
    N <- sum(ns)
    out <- list(NA)
    t1 <- list(NA)
    W <- list(NA)
    for (j in 1:J) {
      W[[j]] <- rep(0, ns[j])
      sum.of.hazards <- rep(0, length(taus[[j]]))
      t1[[j]] <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      sum.of.hazards[1] <- sum(t1[[j]][Rs[[j]][[1]]])
      for (k in 2:length(taus[[j]])) {
        sum.of.hazards[k] <- sum.of.hazards[k-1] - sum(t1[[j]][SDs[[j]][[k-1]]])
      }
      out[[j]] <- sum.of.hazards
    }
    
    # ---------------------------------
    # compute grad
    # ---------------------------------
    Z <- list(NA)
    grad <- list(NA)
    for (j in 1:J) {
      grad[[j]] <- rep(0, ns[j])
      for (k in 1:ns[j]) {
        grad[[j]][k] <- sum(t1[[j]][k]/(out[[j]][Es[[j]] <= k]/Ds[[j]][Es[[j]] <= k])) - dat[[j]]$status[k]
        W[[j]][k] <- sum((t1[[j]][k]*out[[j]][Es[[j]] <= k] - t1[[j]][k]^2)/(out[[j]][Es[[j]] <= k]^2))
      }
      Z[[j]] <- log(t1[[j]]) - grad[[j]]/W[[j]]
    }
    return(list("W" = W, "Z" = Z, "grad" = grad))  
  }
  
  
  ProjCr <- function(beta, r) {
    if (r > 1) {
      temp <- svd(beta)
      beta <- tcrossprod(tcrossprod(temp$u[,1:r, drop=FALSE], diag(temp$d[1:r, drop=FALSE])), temp$v[,1:r, drop=FALSE])
    } else {
      temp <- svd(beta)	
      beta <- tcrossprod(temp$u[,1, drop=FALSE], temp$v[,1, drop=FALSE])*temp$d[1]
    }
    return(beta)
  }
  
  ProjSs <- function(beta, s) {
    if (s < dim(beta)[1]) {
      beta[sort(apply(beta, 1, function(x){sqrt(sum(x^2))}), decreasing=TRUE, index=TRUE)$ix[(s+1):dim(beta)[1]],] <- 0
    }
    return(beta)
  }
  
  
  eval.obj.beta <- function(dat, beta, taus, Es, Rs, Ds, SDs, rho, r, s) {
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    obj <- 0
    for (j in 1:J) {
      sum.of.hazards <- rep(0, length(taus[[j]]))
      t1 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      sum.of.hazards[1] <- sum(t1[Rs[[j]][[1]]])
      obj <- obj - sum(log(t1[Es[[j]][[1]]])) + Ds[[j]][[1]]*log(sum.of.hazards[1])
      for (k in 2:length(taus[[j]])) {
        sum.of.hazards[k] <- sum.of.hazards[k-1] - sum(t1[SDs[[j]][[k-1]]])
        if (!setequal(Es[[j]][[k]],Rs[[j]][[k]])) {
          obj <- obj - sum(log(t1[Es[[j]][[k]]])) + Ds[[j]][k]*log(sum.of.hazards[k])
        }
      }
    }
    return(obj + .25*rho*(sum((beta - ProjCr(beta, r))^2) + sum((beta - ProjSs(beta, s))^2)) + .5*mu*sum(beta^2))
  }

  penalty.method <- function(dat, beta0, taus, Es, Rs, Ds, SDs, rho, r, s) { 
    J <- length(dat)    
    p <- dim(beta0)[1]
    obj.current <- eval.obj.beta(dat, beta0, taus, Es, Rs, Ds, SDs, rho, r, s)      
    for (kk in 1:10) {
      temp <- getGrads(dat, beta0, taus, Ds, Es, Rs, SDs)
      beta.proj <- .5*(ProjCr(beta0, r) + ProjSs(beta0, s))
      beta.temp.1 <- matrix(0, nrow=p, ncol=J)        
      for (j in 1:J) { 
        beta.temp.1[,j] <- updateBeta(xxt=XXts[[j]]/(rho + mu), x=dat[[j]]$x,
        winv = (1/temp$W[[j]]),
        z  = temp$Z[[j]], tildebrho = beta.proj[,j,drop=FALSE]*rho,
        xtw = t(dat[[j]]$x*temp$W[[j]]),
        rhoinv = 1/(rho + mu))
      }
      obj.temp <- eval.obj.beta(dat, beta.temp.1, taus, Es, Rs, Ds, SDs, rho, r, s)
      obj.new <- obj.temp
      beta.new <- beta.temp.1
      if (abs(obj.new - obj.current) < 1e-7*abs(obj.current)) {
        beta0 <- beta.new
        obj.current <- obj.new
        if (!quiet) {
          cat(obj.current, "\n")
        }
        break

      } else {
        beta0 <- beta.new
        obj.current <- obj.new
        if (!quiet) {
         cat(obj.current, "\n")
        }
      }
    }
    return(list("beta" = beta0, "obj.current" = obj.current))
  }

  
  XXts <- list(NA)
  for (j in 1:J) {
    XXts[[j]] <- tcrossprod(dat[[j]]$x)
  }
  betaOut <- array(0, dim=c(p,J,length(svec), length(rvec)))
  beta0 <- matrix(0, nrow=p, ncol=J)
  
  for (s in svec) {
    for (r in rvec) {
      rho <- rho0
      obj.start <- eval.obj.beta(dat, beta0, taus, Es, Rs, Ds, SDs, rho, r, s)
      penMethod <- TRUE
      while (penMethod) {
        update <- try(penalty.method(dat, beta0, taus, Es, Rs, Ds, SDs, rho, r, s))
        if (class(update) == "try-error") {
          rho <- rho0*10
          beta0 <- matrix(0, nrow=p, ncol=J)
        } else {
          beta0 <- update$beta
          obj.new <- update$obj.current
          if (sum((beta0 - ProjCr(beta0, r))^2) + sum((beta0 - ProjSs(beta0, s))^2) < 1e-12) {
            penMethod <- FALSE
          } else {
            rho <- rho*1.2
            if (!silent) {
          	  cat("# ---------------------- ", "\n")
          	  cat("# Through rho = ", rho,  "\n")
          	  cat("# s = ", s, ", r = ", r, "\n")
          	  cat("# ---------------------- ", "\n")
            }
            obj.start <- obj.new
          }
        }
      }
      
      beta0 <- ProjCr(beta0, r)
      beta0 <- ProjSs(beta0, s)
      betaOut[,,which(svec == s), which(rvec == r)] <- beta0
      
    }
  }  
  return(list("beta" = betaOut, "s" = svec, "r" = rvec))
}

