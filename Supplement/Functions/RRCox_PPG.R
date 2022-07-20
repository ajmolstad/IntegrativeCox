library(Rfast)

svd2 <- function(x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE){
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if(!n || !p){
    stop("a dimension is zero")
  }
  La.res <- La.svd(x, nu, nv) 
  res <- list(d = La.res$d)
  if(nu){
    res$u <- La.res$u
  }
  if(nv){
    if(is.complex(x)){
      res$v <- Conj(t(La.res$vt))
    } else {
      res$v <- t(La.res$vt)
    }
  }
  return(res)
}




beta.ProxProxGrad <- function(dat, beta.temp, beta.temp.2, z,  taus, Es, Rs, Ds, lambda, lambda2,  max.iter, tol, inner.quiet){


  J <- length(dat)
  p <- dim(dat[[1]]$x)[2]
  n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
  #W <- matrix(unlist(lapply(dat, function(x){dim(x$x)[1]}))/n, byrow=TRUE, nrow=p, ncol=J)
 # W <- W/max(W)

  eval.obj.beta <- function(dat, beta, beta.svd, taus, Es, Rs, Ds, lambda, lambda2){
    
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    obj <- 0
    for(j in 1:J){
      t1 <- tcrossprod(dat[[j]]$x, t(beta[,j]))
      for(k in 1:length(taus[[j]])){
        obj <- obj - sum(t1[Es[[j]][[k]]]) + 
          Ds[[j]][[k]]*log(sum(exp(t1[Rs[[j]][[k]]])))
      }
    }
    return(obj/n + lambda*sum(beta.svd$d) + #lambda2*sum(W*abs(beta)))
    lambda2*sum(apply(beta, 1, function(x){sqrt(sum(x^2))})))
  }

  eval.loss.beta <- function(dat, beta, taus, Es, Rs, Ds){
    
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    obj <- 0
    for(j in 1:J){
      t1 <- tcrossprod(dat[[j]]$x, t(beta[,j]))
      for(k in 1:length(taus[[j]])){
        obj <- obj - sum(t1[Es[[j]][[k]]]) + 
          Ds[[j]][[k]]*log(sum(exp(t1[Rs[[j]][[k]]])))
      }
    }
    return(obj/n)
  }
    
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs, lambda){

    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){

      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])

      for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
      
    }

    return(gradout/n)
  }

  X.t <- beta.temp
  U.t <- z
  if(is.null(beta.temp.2)){
    Z.t <- beta.temp
  } else {
    Z.t <- beta.temp.2
  }


  rowSoft <- function(x, tau){
    temp <- apply(x, 1, function(x){sqrt(sum(x^2))})
    out <- matrix(0, nrow=dim(x)[1], ncol=dim(x)[2])  
    for(j in which(temp > tau)){
      out[j,] <- x[j,]*(1 - tau/temp[j])
    }
    return(out)
  }

  obj.func <- rep(0, max.iter)
  iterating <- TRUE
  k.iter <- 1
  tau <- 0.75

  # -----------------------------
  # Get initial step size
  # -----------------------------
  eps <- 100
  f.Z0 <- eval.loss.beta(dat, Z.t, taus, Es, Rs, Ds)
  grad.Z0 <- grad.beta(dat, Z.t, taus, Ds, Es, Rs, lambda)
  for(kk in 1:10){
    tildeZ <- Z.t - eps*grad.Z0
    f.Zt <- eval.loss.beta(dat, tildeZ, taus, Es, Rs, Ds)
    if(f.Zt <= f.Z0){
      break
    } else {
      eps <- eps/10
    }
  }
  gamma0 <- 4*(f.Z0 - f.Zt)/sum(grad.Z0^2)

  #gamma0 <- 10
  gamma.t <- gamma0

  while(iterating){

    linesearch <- TRUE
    grad.temp <- grad.beta(dat, Z.t, taus, Ds, Es, Rs, lambda)
    f.Z.t <- eval.loss.beta(dat, Z.t, taus, Es, Rs, Ds)

    while(linesearch){
      X.tp1 <- rowSoft(Z.t - gamma.t*U.t - gamma.t*grad.temp, gamma.t*lambda2)
      #X.tp1 <- pmax(abs(Z.t - gamma.t*U.t - gamma.t*grad.temp) - W*gamma.t*lambda2, 0)*sign(Z.t - gamma.t*U.t - gamma.t*grad.temp)
      #eo <- svd(Z.t - gamma.t*U.t - gamma.t*grad.temp)
      #X.tp1 <- tcrossprod(eo$u*(tcrossprod(rep(1,dim(eo$u)[1]), pmax(eo$d - gamma.t*lambda, 0))),eo$v)
      f.X.tp1 <- eval.loss.beta(dat, X.tp1, taus, Es, Rs, Ds)
      Q.X.tp1 <- f.Z.t + sum(grad.temp*(X.tp1 - Z.t)) + (sum((X.tp1 - Z.t)^2)/(2*gamma.t))
      if(f.X.tp1 > Q.X.tp1){
        gamma.t <- tau*gamma.t
      } else {
        linesearch <- FALSE
      }
    }

    #Z.tp1 <- rowSoft(X.tp1 + gamma.t*U.t, gamma.t*lambda2)
    eo <- svd2(X.tp1 + gamma.t*U.t)
    Z.tp1 <- tcrossprod(eo$u*(tcrossprod(rep(1,dim(eo$u)[1]), pmax(eo$d - gamma.t*lambda, 0))),eo$v)
    U.tp1 <- U.t + (X.tp1 - Z.tp1)/gamma.t

    # ----------------------------------------
    # Get next step size and update iterates
    # ----------------------------------------
    # delta.t <- Q.X.tp1 - f.X.tp1
    # gamma.t <- sqrt(gamma.t^2 + (gamma.t*delta.t)/(4*p*J*lambda))

    # if(k.iter %% 10 == 0){
    #   gamma.t <- gamma0
    # }

    Z.t <- Z.tp1
    X.t <- X.tp1
    U.t <- U.tp1


    # ----------------------------------------
    # Check convergence
    # ----------------------------------------
    if(!inner.quiet){
      cat(obj.func[k.iter] <- eval.obj.beta(dat, X.tp1, svd2(X.tp1), taus, Es, Rs, Ds, lambda, lambda2), "\n")
    } else {
      obj.func[k.iter] <- eval.obj.beta(dat, X.tp1, svd2(X.tp1), taus, Es, Rs, Ds, lambda, lambda2)
    }
    if(k.iter > 3){
      if(abs(obj.func[k.iter-1] - obj.func[k.iter]) < tol & abs(obj.func[k.iter-2] - obj.func[k.iter-1]) < tol & 
      abs(obj.func[k.iter-3] - obj.func[k.iter-2]) < tol){
        iterating <- FALSE
      } 
    }

    if(k.iter > max.iter){
      iterating <- FALSE
    }

    k.iter <- k.iter + 1
    # cat(gamma.t, "\n")
    # cat("# --------------------- ")

  }


  # beta.koh <- rowSoft(z, alpha*lambda2)
  # b.Out <- X.tp1
  # X.tp1[which(apply(Xtp1, 1, function(x){sqrt(sum(x^2))}) < 1e-6),] <- 0
  # getVals <- svd(Z.tp1)
  # temp <- (1*(X.tp1!=0))*Z.tp1
  # eo <- svd(temp)
  # if(any(abs(getVals$d) > 1e-8)){
  #   r <- max(which(abs(getVals$d) > 1e-8))
  #   beta.out <- eo$u[,1:r,drop=FALSE]%*%diag(eo$d[1:r])%*%t(eo$v[,1:r, drop=FALSE])
  # } else {
  #   beta.out <- matrix(0, nrow=p, ncol=J)
  # }
  return(list("beta" = (1*(X.tp1!=0))*Z.tp1, "z" = U.tp1, "beta.2" = Z.tp1))
  #return(list("beta" = Z.tp1, "z" = U.tp1, "beta.2" = X.tp1))
  # checkMat <- Z.tp1
  # for(j in 1:p){
  #   checkMat[j,] <- Z.tp1[j,]/sqrt(sum(Z.tp1[j,]^2))
  # }
  # grad.beta(dat, X.tp1, taus, Ds, Es, Rs, lambda) + lambda*svd(X.tp1)$u%*%t(svd(X.tp1)$v) + lambda2*checkMat
 
}


beta.AccProxGrad.SVD <- function(dat, beta.temp, taus, Es, Rs, Ds, lambda, max.iter, tol, inner.quiet){

  J <- length(dat)
  p <- dim(dat[[1]]$x)[2]
  n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
  
  eval.obj.beta <- function(dat, beta, beta.svd, taus, Es, Rs, Ds, lambda){
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    obj <- 0
    for(j in 1:J){
      t1 <- tcrossprod(dat[[j]]$x, t(beta[,j]))
      for(k in 1:length(taus[[j]])){
        obj <- obj - sum(t1[Es[[j]][[k]]]) + 
          Ds[[j]][[k]]*log(sum(exp(t1[Rs[[j]][[k]]])))
      }
    }
    return(obj/n + lambda*sum(beta.svd$d))
  }
    
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs){

    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){

      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])

     for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
    }

    return(gradout/n)
  }

  beta.km1 <- beta.temp
  beta.km2 <- beta.km1

  alphakm1 <- 1
  alphakm2 <- 1
  L0 <- n/max(sapply(1:J, function(x){max(svd(dat[[x]]$x)$d)^2}))
  obj.func <- rep(0, max.iter)
  iterating <- TRUE
  k.iter <- 1
  pll <- rep(0, max.iter)

  while(iterating){

    # --------------------
    # Update
    # --------------------
    L <- L0*200
    Atemp <- beta.km1 + ((alphakm2 - 1)/alphakm1)*(beta.km1 - beta.km2)
    grad <- grad.beta(dat, Atemp, taus, Ds, Es, Rs)
    lik <- eval.obj.beta(dat, Atemp, svd2(Atemp), taus, Es, Rs, Ds, lambda = 0)
    linesearch <- TRUE
    
    # ----------------------------
    # Proximal update
    # ----------------------------
    while(linesearch){
      
      beta.svd <- svd2(Atemp - L*grad)
      beta.temp <- beta.svd$u%*%diag(pmax(beta.svd$d - L*lambda, 0))%*%t(beta.svd$v)
      beta.svd$d <- pmax(beta.svd$d - L*lambda, 0)
      templik <- eval.obj.beta(dat, beta.temp, beta.svd, taus, Es, Rs, Ds, lambda = 0) 

     if(L == L0){

        linesearch <- FALSE
        pll[k.iter] <- eval.obj.beta(dat, beta.temp, beta.svd, taus, Es, Rs, Ds, lambda) 
        beta.km2 <- beta.km1
        beta.k <- beta.temp

     } else {
        
        if(templik < (lik + sum(diag(crossprod(grad, beta.temp - Atemp))) + 1/(2*L)*sum((beta.temp - Atemp)^2))){
            
            linesearch <- FALSE
            beta.k <- beta.temp
            beta.km2 <- beta.km1
            pll[k.iter] <- eval.obj.beta(dat, beta.k, beta.svd, taus, Es, Rs, Ds, lambda) 
            
        } else {
          L <- max(L/2, L0)
        }
     }
      
   }

    
    beta.km1 <- beta.k
    alphakm2 <- alphakm1
    alphakm1 <- (1 + sqrt(1 + 4*alphakm1^2))/2
    if(!inner.quiet){
      cat("# ------------- ", "\n")
      cat(pll[k.iter], "\n")
    }

    if(k.iter> 3){
      if(abs(pll[k.iter] - pll[k.iter-1]) < tol && abs(pll[k.iter-1] - pll[k.iter-2]) < tol && abs(pll[k.iter-2] - pll[k.iter-3]) < tol){
        iterating <- FALSE
      }
    }

    k.iter <- k.iter + 1
    if(k.iter > max.iter){
      iterating <- FALSE
    }

  }



  return(list("beta" = beta.k))

 }

beta.AccProxGrad.L2 <- function(dat, beta.temp, taus, Es, Rs, Ds, lambda2, max.iter, tol, inner.quiet){

  J <- length(dat)
  p <- dim(dat[[1]]$x)[2]
  alpha <- 1
  n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
  
  eval.obj.beta <- function(dat, beta, taus, Es, Rs, Ds, lambda2){
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    obj <- 0
    for(j in 1:J){
      t1 <- tcrossprod(dat[[j]]$x, t(beta[,j]))
      for(k in 1:length(taus[[j]])){
        obj <- obj - sum(t1[Es[[j]][[k]]]) + 
          Ds[[j]][[k]]*log(sum(exp(t1[Rs[[j]][[k]]])))
      }
    }
    return(obj/n + lambda2*sum(apply(beta, 1, function(x){sqrt(sum(x^2))})))
  }
    
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs){

    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){

      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])

     for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
    }

    return(gradout/n)
  }


  rowSoft <- function(x, tau){
    temp <- apply(x, 1, function(x){sqrt(sum(x^2))})
    out <- matrix(0, nrow=dim(x)[1], ncol=dim(x)[2])  
    for(j in which(temp > tau)){
      out[j,] <- x[j,]*(1 - tau/temp[j])
    }
    return(out)
  }


  beta.km1 <- beta.temp
  beta.km2 <- beta.km1

  alphakm1 <- 1
  alphakm2 <- 1
  L0 <- n/max(sapply(1:J, function(x){max(svd2(dat[[x]]$x)$d)^2}))
  obj.func <- rep(0, max.iter)
  iterating <- TRUE
  k.iter <- 1
  pll <- rep(0, max.iter)

  while(iterating){

    # --------------------
    # Update
    # --------------------
    L <- L0*200
    Atemp <- beta.km1 + ((alphakm2 - 1)/alphakm1)*(beta.km1 - beta.km2)
    grad <- grad.beta(dat, Atemp, taus, Ds, Es, Rs)
    lik <- eval.obj.beta(dat, Atemp, taus, Es, Rs, Ds, lambda2 = 0)
    linesearch <- TRUE
    
    # ----------------------------
    # Proximal update
    # ----------------------------
    while(linesearch){
      
      beta.temp <- rowSoft(Atemp - L*grad, L*lambda2)
      templik <- eval.obj.beta(dat, beta.temp, taus, Es, Rs, Ds, lambda2 = 0) 

     if(L == L0){

        linesearch <- FALSE
        pll[k.iter] <- eval.obj.beta(dat, beta.temp, taus, Es, Rs, Ds, lambda2) 
        beta.km2 <- beta.km1
        beta.k <- beta.temp

     } else {
        
        if(templik < (lik + sum(diag(crossprod(grad, beta.temp - Atemp))) + 1/(2*L)*sum((beta.temp - Atemp)^2))){
            
            linesearch <- FALSE
            beta.k <- beta.temp
            beta.km2 <- beta.km1
            pll[k.iter] <- eval.obj.beta(dat, beta.k, taus, Es, Rs, Ds, lambda2) 
            
        } else {
          L <- max(L/2, L0)
        }
     }
      
   }

    
    beta.km1 <- beta.k
    alphakm2 <- alphakm1
    alphakm1 <- (1 + sqrt(1 + 4*alphakm1^2))/2
    if(!inner.quiet){
      cat("# ------------- ", "\n")
      cat(pll[k.iter], "\n")
    }

    if(k.iter> 3){
      if(abs(pll[k.iter] - pll[k.iter-1]) < tol && abs(pll[k.iter-1] - pll[k.iter-2]) < tol && abs(pll[k.iter-2] - pll[k.iter-3]) < tol){
        iterating <- FALSE
      }
    }

    k.iter <- k.iter + 1
    if(k.iter > max.iter){
      iterating <- FALSE
    }

  }



  return(list("beta" = beta.k))

 }


RRCox.coef <- function(fit, lambda = NULL){
  
  if(is.null(lambda)){
    warning("Need to input tuning parameter lambda!", "\n")
  }
  beta <- fit$beta[which(fit$lambda.vec == lambda),,]
  if(fit$standardize){
    message("Note that output coefficients are on standardized predictor scale.")
  }
  return(beta)
  
}

RRCox.predict <- function(fit, newX = NULL, lambda = NULL, which.pop = NULL){
  
  if(is.null(lambda)){
    warning("Need to input tuning parameter lambda!", "\n")
  }
  if(is.null(newX)){
    warning("Need to input newX values!", "\n")
  }
  if(is.null(which.pop)){
    warning("Need to specify which population to predict!", "\n")
  }
  
  if(fit$standardize){
    x <- (newX - rep(1, dim(newX)[1])%*%t(fit$X.mean[[which.pop]]))/(rep(1, dim(newX)[1])%*%t(fit$X.sd[[which.pop]]))
  } else {
    x <- newX
  }
  preds <- x%*%fit$beta[which(fit$lambda.vec==lambda),,which.pop]
  return(preds)
  
}


# --------------------------------------------
# Reduced rank cox model CV function
# --------------------------------------------
RRCox.cv <- function(dat, 
                     nlambda1 = 100, nlambda2 = 0,
                     delta1 = .2, delta2 = .2, 
                     nfolds = NULL, 
                     standardize.x = FALSE,
                     max.iter = 500,
                     inner.quiet = TRUE, 
                     quiet = TRUE,
                     inner.tol = 1e-8){
  
  if(nlambda2 > 0 & nlambda1 > 0){
        out <- RRCox.cv.lam1lam2(dat = dat, nlambda1 = nlambda1, nlambda2 = nlambda2,
                                             delta1 = delta1, delta2 = delta2, 
                                             nfolds = nfolds, 
                                             standardize.x = standardize.x,
                                             max.iter = max.iter,
                                             inner.quiet = inner.quiet, 
                                             quiet = quiet,
                                             inner.tol = inner.tol)
  } else {
    if(nlambda2 == 0){
      out <- RRCox.cv.lam1(dat = dat, nlambda1 = nlambda1,
                              delta1 = delta1,
                              nfolds = nfolds, 
                              standardize.x = standardize.x,
                              max.iter = max.iter,
                              inner.quiet = inner.quiet, 
                              quiet = quiet,
                              inner.tol = inner.tol)
      } else {
        out <- RRCox.cv.lam2(dat = dat, nlambda2 = nlambda2,
                              delta2 = delta2,
                              nfolds = nfolds, 
                              standardize.x = standardize.x,
                              max.iter = max.iter,
                              inner.quiet = inner.quiet, 
                              quiet = quiet,
                              inner.tol = inner.tol)
      }
  }
  
  return(out)
}
  
  
  

# -- if nlambda1 > 1 and nlambda2 > 0
RRCox.cv.lam1lam2 <- function(dat, 
                    nlambda1, nlambda2,
                    delta1, delta2, 
                    nfolds = NULL, 
                    standardize.x = FALSE,
                    max.iter = 500,
                    inner.quiet = TRUE, 
                    quiet = TRUE,
                    inner.tol = 1e-8){
  
  
  # ------------------------------------------
  # Preliminaries
  # ------------------------------------------
  J <- length(dat)
  q <- dim(dat[[1]]$z)[2]

  taus <- list(NA)
  Es <- list(NA)
  Rs <- list(NA)
  Ds <- list(NA)
  for(jj in 1:J){
    taus[[jj]] <- dat[[jj]]$time[which(dat[[jj]]$status==1)]
    Es[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time == x)})
    Ds[[jj]] <- sapply(taus[[jj]], function(x){length(which(dat[[jj]]$time == x))})
    Rs[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time >= x)})
  }

  n <- sum(unlist(lapply(dat, function(x){dim(x$X)[1]})))

  # if(!is.null(nfolds)){
  #   cv.index <- list(NA)
  #   for(j in 1:J){
  #     fold1 <- sample(rep(1:nfolds, length=dim(dat[[j]]$X)[1]))
  #     cv.index[[j]] <- split(1:dim(dat[[j]]$X)[1], fold1)
  #   }
  # }

  if(!is.null(nfolds)){
    cv.index <- list(NA)
    for(j in 1:J){
      fold1 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==1)))
      uncens <- which(dat[[j]]$status==1)
      cv.index[[j]] <- split(uncens, fold1)

      fold0 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==0)))
      cens <- which(dat[[j]]$status==0)
      temp0 <- split(cens, fold0)
      for(jj in 1:nfolds){
        cv.index[[j]][[jj]] <- c(cv.index[[j]][[jj]], temp0[[jj]])
      }
    }
  }
  
  
  # ------------------------------------------------------
  # Remove genes with no variability in at least one fold
  # ------------------------------------------------------
  if(!is.null(nfolds)){
    rm.genes <- NULL
    for(kk in 1:nfolds){
      for(j in 1:J){
        temp <- which(apply(dat[[j]]$X[cv.index[[j]][[kk]],], 2, sd) == 0)
        if(length(temp)> 0){
          rm.genes <- c(rm.genes, temp)
        }
      }
    }
  } else {
    rm.genes <- NULL
    for(j in 1:J){
      temp <- which(apply(dat[[j]]$X, 2, sd) == 0)
      if(length(temp)> 0){
        rm.genes <- c(rm.genes, temp)
      }
    }
  }
  
  if(length(rm.genes) > 0){
    for(j in 1:J){
      dat[[j]]$X <- dat[[j]]$X[,-rm.genes]
    }
  }
  p <- dim(dat[[1]]$X)[2]
  
  # ------------------------------------------------------
  # Standardize.x if necessary
  # ------------------------------------------------------
  if(standardize.x){
    # --- first check if need to 
    for(jj in 1:J){
      dat[[jj]]$x <- (dat[[jj]]$X - rep(1, dim(dat[[jj]]$X)[1])%*%t(colMeans(dat[[jj]]$X)))/(rep(1, dim(dat[[jj]]$X)[1])%*%t(apply(dat[[jj]]$X, 2, sd)))
    }
  } else {
    for(jj in 1:J){
      dat[[jj]]$x <- dat[[jj]]$X
    }
  }
  
  
  # ------------------------------------------------------
  # Get tuning candidate tuning parameter range 
  # ------------------------------------------------------
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs, lambda){

    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){

      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])

     for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
    }

    return(gradout/n)
  }

  ZetaList <- list(NA)
  beta.init <- matrix(0, nrow=p, ncol=J)
  for(j in 1:J){
    ZetaList[[j]] <- rep(0, dim(dat[[j]]$x)[1])
  }
  tempGrad <- grad.beta(dat, beta.init, taus, Ds, Es, Rs, lambda)
  
  lam.max <- max(apply(tempGrad,1, function(x){sqrt(sum(x^2))}))
  lam.min <- delta2*lam.max
  lambda2.vec <- 10^seq(log10(lam.max), log10(lam.min), length=nlambda2)

  # ------------------------------------------------------
  # Output matrices 
  # ------------------------------------------------------
  beta.out <- array(0, dim=c(nlambda1+1, nlambda2, p, J))
  
  # ------------------------------------------------
  # Generate necessary preliminary quantities
  # ------------------------------------------------
  beta.temp <- beta.init
  z.init <- matrix(0, nrow=p, ncol=J)
  z.temp <- z.init

  for(kk in 1:length(lambda2.vec)){
    
    temp <- beta.AccProxGrad.L2(dat = dat, beta.temp = beta.temp, taus = taus, Es = Es, 
                          Rs = Rs, Ds = Ds, lambda2 = lambda2.vec[kk], 
                          max.iter = max.iter, tol = inner.tol,
                          inner.quiet = inner.quiet)
    beta.temp <- temp$beta
    beta.out[1,kk,,] <- beta.temp
    cat(kk, "\n")
    cat("# -------------- ", "\n")
  }

  lambda1.vec <- rep(0, nlambda1+1)
  tempGrad <- grad.beta(dat, matrix(0, nrow=p, ncol=J), taus, Ds, Es, Rs, lambda)
  lam1.max <- max(svd(tempGrad)$d)
  lam1.min <- delta1*lam1.max
  for(ll in 1:nlambda1){
    lambda1.vec[ll] <- (lam1.max^((nlambda1-ll)/(nlambda1-1)))*(lam1.min^((ll-1)/(nlambda1-1)))
  }
  lambda1.vec <- rev(lambda1.vec)


  for(kk in 1:length(lambda2.vec)){

      beta.temp <- beta.out[1,kk,,]
      z.temp <- matrix(0, nrow=p, ncol=J)
      beta.temp.2 <- beta.out[1,kk,,]

        for(zz in 2:length(lambda1.vec)){ 
          temp <- beta.ProxProxGrad(dat = dat, beta.temp = beta.temp, beta.temp.2 = beta.temp.2, z = z.temp,
                                   taus = taus, Es = Es, Rs = Rs, Ds = Ds, 
                                   lambda = lambda1.vec[zz], lambda2.vec[kk], 
                                   max.iter = max.iter, tol = inner.tol, 
                                   inner.quiet = inner.quiet)

          beta.temp <- temp$beta
          z.temp <- temp$z
          beta.temp.2 <- temp$beta.2
          if(zz == 2){
            z.init <- z.temp
          }
          beta.out[zz,kk,,] <- beta.temp
          if(!quiet){
            cat("# ------------------------------------------------------", "\n")
            cat("Through tuning parameter", kk, "/", zz, " of ", length(lambda2.vec),"/", length(lambda1.vec),  "\n")
            cat(sum(rowSums(beta.temp!=0) > 0),":", sum(beta.temp^2), "\n")
          }

          if(sum(beta.temp==0) == dim(beta.temp)[1]*dim(beta.temp)[2]){
            break
          }

        }

    }
      
  if(!is.null(nfolds)){
    
    cat("Into CV", "\n")
    # -----------------------------------------
    # save metric for cross-validation 
    # -----------------------------------------
    PartLik <- array(NA, dim=c(length(lambda1.vec), length(lambda2.vec), nfolds, J))
    linPredTest <- list(NA)
    for(j in 1:J){
      linPredTest[[j]] <- array(0, dim = c(dim(dat[[j]]$X)[1], length(lambda1.vec), length(lambda2.vec)))
    }
    linPredPartLik <- array(0, dim=c(length(lambda1.vec), length(lambda2.vec), J))
      
    for(kk in 1:nfolds){
      
      tausTrain <- list(NA)
      EsTrain <- list(NA)
      RsTrain <- list(NA)
      DsTrain <- list(NA)
      taus <- NULL
      Es <- NULL
      Rs <- NULL
      Ds <- NULL
      
      # ------------------------------------------------
      # Standardize and creating training/testing data
      # ------------------------------------------------
      for(jj in 1:J){
        
        tempTrain <- dat[[jj]]$X[-cv.index[[jj]][[kk]],]
        tempTest <- dat[[jj]]$X[cv.index[[jj]][[kk]],]
        if(!standardize.x){
          dat[[jj]]$x <- tempTrain
          dat[[jj]]$xTest <- tempTest
        } else {
          dat[[jj]]$x <- (tempTrain - rep(1, dim(tempTrain)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTrain)[1])%*%t(apply(tempTrain, 2, sd)))
          dat[[jj]]$xTest <- (tempTest - rep(1, dim(tempTest)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTest)[1])%*%t(apply(tempTrain, 2, sd)))
        }
        dat[[jj]]$timeTrain <- dat[[jj]]$time[-cv.index[[jj]][[kk]]]
        dat[[jj]]$timeTest <- dat[[jj]]$time[cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTrain <- dat[[jj]]$status[-cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTest <-  dat[[jj]]$status[cv.index[[jj]][[kk]]]
        tausTrain[[jj]] <- dat[[jj]]$timeTrain[which(dat[[jj]]$statusTrain==1)]
        EsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain == x)})
        DsTrain[[jj]] <- sapply(tausTrain[[jj]],function(x){length(which(dat[[jj]]$timeTrain == x))})
        RsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain >= x)})
        
      }

      # ------------------------------------------------
      # Generate necessary preliminary quantities
      # ------------------------------------------------
      n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
      beta.temp.2 <- NULL
      z.init <- matrix(0, nrow=p, ncol=J)
      z.temp <- z.init
      beta.init.store <- array(0, dim=c(p, J, length(lambda2.vec)))
      for(lam2ind in 1:length(lambda2.vec)){

        temp0 <- beta.AccProxGrad.L2(dat = dat, beta.temp = beta.temp, taus = tausTrain, Es = EsTrain, 
                              Rs = RsTrain, Ds = DsTrain, lambda2 = lambda2.vec[lam2ind], 
                              max.iter = max.iter, tol = inner.tol,
                              inner.quiet = inner.quiet)
        beta.temp <- temp0$beta
        beta.init.store[,,lam2ind] <- temp0$beta

        for(j in 1:J){
          dev1 <- coxnet.deviance(pred = dat[[j]]$x%*%temp0$beta[,j], y = Surv(dat[[j]]$timeTrain, dat[[j]]$statusTrain))
          dev2 <- coxnet.deviance(pred = dat[[j]]$X%*%temp0$beta[,j], y = Surv(dat[[j]]$time, dat[[j]]$status))
          PartLik[1, lam2ind,  kk, j] <- dev2 - dev1
          linPredTest[[j]][cv.index[[j]][[kk]],1,lam2ind] <- dat[[j]]$xTest%*%temp0$beta[,j]
        } 

        if(!inner.quiet){
          cat("Through ",1, "/", lam2ind, " of ", length(lambda1.vec), "/", length(lambda2.vec), "in fold", kk, "\n")
        }

      }

      for(lam2ind in 1:length(lambda2.vec)){

        beta.temp <- beta.init.store[,,lam2ind]
        z.temp <- matrix(0, nrow=p, ncol=J)
        beta.temp.2 <- beta.init.store[,,lam2ind]

        for(lam1ind in 2:length(lambda1.vec)){ 

          temp0 <- beta.ProxProxGrad(dat = dat, beta.temp = beta.temp, beta.temp.2 = beta.temp.2, z = z.temp,
                                   taus = tausTrain, Es = EsTrain, Rs = RsTrain, Ds = DsTrain, 
                                   lambda = lambda1.vec[lam1ind], lambda2.vec[lam2ind], 
                                   max.iter = max.iter, tol = inner.tol, 
                                   inner.quiet = inner.quiet)
          beta.temp <- temp0$beta
          z.temp <- temp0$z
          beta.temp.2 <- temp0$beta.2

          for(j in 1:J){
            dev1 <- coxnet.deviance(pred = dat[[j]]$x%*%temp0$beta[,j], y = Surv(dat[[j]]$timeTrain, dat[[j]]$statusTrain))
            dev2 <- coxnet.deviance(pred = dat[[j]]$X%*%temp0$beta[,j], y = Surv(dat[[j]]$time, dat[[j]]$status))
            PartLik[lam1ind, lam2ind,  kk, j] <- dev2 - dev1
            linPredTest[[j]][cv.index[[j]][[kk]],lam1ind, lam2ind] <- dat[[j]]$xTest%*%temp0$beta[,j]
          } 

          if(!inner.quiet){
            cat("Through ", lam1ind, "/", lam2ind, " of ", length(lambda1.vec), "/", length(lambda2.vec), "in fold", kk, "\n")
          }

          if(sum(beta.temp==0) == dim(beta.temp)[1]*dim(beta.temp)[2]){
            break
          }

        }
      }
    }



    for(ll in 1:length(lambda1.vec)){
      for(lll in 1:length(lambda2.vec)){
        for(jjj in 1:J){
          linPredPartLik[ll,lll,jjj] <- coxnet.deviance(y=Surv(dat[[jjj]]$time, dat[[jjj]]$status), pred=linPredTest[[jjj]][,ll,lll])
        }
      }
    }
  
    
  } else {
    PartLik <- NULL
    linPredPartLik <- NULL
    linPredTest <- NULL
  }
  
  return(list("PartLik" = PartLik, 
              "linPredTest" = linPredTest,
              "linPredPartLik" = linPredPartLik,
              "beta" = beta.out, 
              "rm.genes" = rm.genes,
              "X.mean" = lapply(dat, function(x){colMeans(x$X)}),
              "X.sd" = lapply(dat, function(x){apply(x$X, 2, sd)}),
              "lambda.vec" = lambda1.vec, 
              "lambda2.vec" = lambda2.vec, 
              "standardize" = standardize.x))
  
}


# -- if nlambda1 > 1 and nlambda2 == 0 
RRCox.cv.lam1 <- function(dat,nlambda1,
                              delta1, 
                              nfolds = NULL, 
                              standardize.x = FALSE,
                              max.iter = 500,
                              inner.quiet = TRUE, 
                              quiet = TRUE,
                              inner.tol = 1e-8){
  
  
  # ------------------------------------------
  # Preliminaries
  # ------------------------------------------
  J <- length(dat)
  q <- dim(dat[[1]]$z)[2]
  n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
  
  taus <- list(NA)
  Es <- list(NA)
  Rs <- list(NA)
  Ds <- list(NA)
  for(jj in 1:J){
    taus[[jj]] <- dat[[jj]]$time[which(dat[[jj]]$status==1)]
    Es[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time == x)})
    Ds[[jj]] <- sapply(taus[[jj]], function(x){length(which(dat[[jj]]$time == x))})
    Rs[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time >= x)})
  }
  

  if(!is.null(nfolds)){
    cv.index <- list(NA)
    for(j in 1:J){
      fold1 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==1)))
      uncens <- which(dat[[j]]$status==1)
      cv.index[[j]] <- split(uncens, fold1)
      
      fold0 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==0)))
      cens <- which(dat[[j]]$status==0)
      temp0 <- split(cens, fold0)
      for(jj in 1:nfolds){
        cv.index[[j]][[jj]] <- c(cv.index[[j]][[jj]], temp0[[jj]])
      }
    }
  }
  
  
  # ------------------------------------------------------
  # Remove genes with no variability in at least one fold
  # ------------------------------------------------------
  if(!is.null(nfolds)){
    rm.genes <- NULL
    for(kk in 1:nfolds){
      for(j in 1:J){
        temp <- which(apply(dat[[j]]$X[cv.index[[j]][[kk]],], 2, sd) == 0)
        if(length(temp)> 0){
          rm.genes <- c(rm.genes, temp)
        }
      }
    }
  } else {
    rm.genes <- NULL
    for(j in 1:J){
      temp <- which(apply(dat[[j]]$X, 2, sd) == 0)
      if(length(temp)> 0){
        rm.genes <- c(rm.genes, temp)
      }
    }
  }
  
  if(length(rm.genes) > 0){
    for(j in 1:J){
      dat[[j]]$X <- dat[[j]]$X[,-rm.genes]
    }
  }
  p <- dim(dat[[1]]$X)[2]
  
  # ------------------------------------------------------
  # Standardize.x if necessary
  # ------------------------------------------------------
  if(standardize.x){
    # --- first check if need to 
    for(jj in 1:J){
      dat[[jj]]$x <- (dat[[jj]]$X - rep(1, dim(dat[[jj]]$X)[1])%*%t(colMeans(dat[[jj]]$X)))/(rep(1, dim(dat[[jj]]$X)[1])%*%t(apply(dat[[jj]]$X, 2, sd)))
    }
  } else {
    for(jj in 1:J){
      dat[[jj]]$x <- dat[[jj]]$X
    }
  }
  
  
  # ------------------------------------------------------
  # Get tuning candidate tuning parameter range 
  # ------------------------------------------------------
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs, lambda){
    
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){
      
      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])
      
      for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
    }
    
    return(gradout/n)
  }
  
  ZetaList <- list(NA)
  beta.init <- matrix(0, nrow=p, ncol=J)
  for(j in 1:J){
    ZetaList[[j]] <- rep(0, dim(dat[[j]]$x)[1])
  }
  tempGrad <- grad.beta(dat, beta.init, taus, Ds, Es, Rs, lambda)

  lam.max <- max(svd2(tempGrad)$d)
  lam.min <- delta1*lam.max
  lambda.vec <- 10^seq(log10(lam.max), log10(lam.min), length=nlambda1)
    
  
  # ------------------------------------------------------
  # Output matrices 
  # ------------------------------------------------------
  beta.out <- array(0, dim=c(length(lambda.vec), 1, p, J))
  
  # ------------------------------------------------
  # Generate necessary preliminary quantities
  # ------------------------------------------------
  beta.temp <- beta.init

  for(kk in 1:length(lambda.vec)){
    
    temp <- beta.AccProxGrad.SVD(dat = dat, beta.temp = beta.temp, taus = taus, Es = Es, 
                                Rs = Rs, Ds = Ds, lambda = lambda.vec[kk], 
                                max.iter = max.iter, tol = inner.tol, 
                                inner.quiet =   inner.quiet)
    beta.temp <- temp$beta
    beta.init <- beta.temp
    
    if(!quiet){
      cat("# ------------------------------------------------------", "\n")
      cat("Through tuning parameter", kk, " of ", length(lambda.vec),  "\n")
      cat(sum(rowSums(beta.temp!=0) > 0),":", sum(beta.temp^2), "\n")
    }
    
    beta.out[kk,1,,] <- beta.temp
    
  }

  
  if(!is.null(nfolds)){
    
    # -----------------------------------------
    # save metric for cross-validation 
    # -----------------------------------------
    PartLik <- array(NA, dim=c(length(lambda.vec), 1, nfolds, J))
    linPredTest <- list(NA)
    for(j in 1:J){
      linPredTest[[j]] <- array(0, dim = c(dim(dat[[j]]$X)[1], length(lambda.vec), 1))
    }
    linPredPartLik <- array(0, dim=c(length(lambda.vec), 1, J))
    
    for(kk in 1:nfolds){
    
      tausTrain <- list(NA)
      EsTrain <- list(NA)
      RsTrain <- list(NA)
      DsTrain <- list(NA)
      taus <- NULL
      Es <- NULL
      Rs <- NULL
      Ds <- NULL
      
      # ------------------------------------------------
      # Standardize and creating training/testing data
      # ------------------------------------------------
      for(jj in 1:J){
        
        tempTrain <- dat[[jj]]$X[-cv.index[[jj]][[kk]],]
        tempTest <- dat[[jj]]$X[cv.index[[jj]][[kk]],]
        if(!standardize.x){
          dat[[jj]]$x <- tempTrain
          dat[[jj]]$xTest <- tempTest
        } else {
          dat[[jj]]$x <- (tempTrain - rep(1, dim(tempTrain)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTrain)[1])%*%t(apply(tempTrain, 2, sd)))
          dat[[jj]]$xTest <- (tempTest - rep(1, dim(tempTest)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTest)[1])%*%t(apply(tempTrain, 2, sd)))
        }
        dat[[jj]]$timeTrain <- dat[[jj]]$time[-cv.index[[jj]][[kk]]]
        dat[[jj]]$timeTest <- dat[[jj]]$time[cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTrain <- dat[[jj]]$status[-cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTest <-  dat[[jj]]$status[cv.index[[jj]][[kk]]]
        tausTrain[[jj]] <- dat[[jj]]$timeTrain[which(dat[[jj]]$statusTrain==1)]
        EsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain == x)})
        DsTrain[[jj]] <- sapply(tausTrain[[jj]],function(x){length(which(dat[[jj]]$timeTrain == x))})
        RsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain >= x)})
        
      }
      
      # -----------------------------------------------
      # Get tuning candidate tuning parameter range 
      # -----------------------------------------------
      beta.temp <- matrix(0, nrow=p, ncol=J)

      for(ll in 1:length(lambda.vec)){
        
        temp0 <- beta.AccProxGrad.SVD(dat = dat, beta.temp = beta.temp, taus = tausTrain, Es = EsTrain, 
                                     Rs = RsTrain, Ds = DsTrain, lambda = lambda.vec[ll], 
                                     max.iter = max.iter, tol = inner.tol, 
                                     inner.quiet = TRUE)
        beta.temp <- temp0$beta
        if(!inner.quiet){
          cat("Through ", ll, " of ", length(lambda.vec), " in fold ", kk, " of ", nfolds, "\n")
        }

        for(j in 1:J){
          dev1 <- coxnet.deviance(pred = dat[[j]]$x%*%temp0$beta[,j], y = Surv(dat[[j]]$timeTrain, dat[[j]]$statusTrain))
          dev2 <- coxnet.deviance(pred = dat[[j]]$X%*%temp0$beta[,j], y = Surv(dat[[j]]$time, dat[[j]]$status))
          PartLik[ll,1,  kk, j] <- dev2 - dev1
          linPredTest[[j]][cv.index[[j]][[kk]],ll,1] <- dat[[j]]$xTest%*%temp0$beta[,j]
        }
      }
      
      if(!quiet){
        cat("# ------------------------- ", "\n")
        cat("Through fold ", kk, " of ", nfolds, "\n")
      }
    }
    
    for(ll in 1:length(lambda.vec)){
      for(jjj in 1:J){
        linPredPartLik[ll,1,jjj] <- coxnet.deviance(y=Surv(dat[[jjj]]$time, dat[[jjj]]$status), pred=linPredTest[[jjj]][,ll,1])
      }
    }
    
    
  } else {
    PartLik <- NULL
    linPredPartLik <- NULL
    linPredTest <- NULL
  }
  
  
  
  
  
  return(list("PartLik" = PartLik, 
              "linPredTest" = linPredTest,
              "linPredPartLik" = linPredPartLik,
              "beta" = beta.out, 
              "rm.genes" = rm.genes,
              "X.mean" = lapply(dat, function(x){colMeans(x$X)}),
              "X.sd" = lapply(dat, function(x){apply(x$X, 2, sd)}),
              "lambda.vec" = lambda.vec, 
              "lambda2.vec" = 0, 
              "standardize" = standardize.x))
  
}












# -- if nlambda1 > 1 and nlambda2 == 0 
RRCox.cv.lam2 <- function(dat,nlambda2,
                              delta2, 
                              nfolds = NULL, 
                              standardize.x = FALSE,
                              max.iter = 500,
                              inner.quiet = TRUE, 
                              quiet = TRUE,
                              inner.tol = 1e-8){
  
  
  # ------------------------------------------
  # Preliminaries
  # ------------------------------------------
  J <- length(dat)
  q <- dim(dat[[1]]$z)[2]
  n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
  
  taus <- list(NA)
  Es <- list(NA)
  Rs <- list(NA)
  Ds <- list(NA)
  for(jj in 1:J){
    taus[[jj]] <- dat[[jj]]$time[which(dat[[jj]]$status==1)]
    Es[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time == x)})
    Ds[[jj]] <- sapply(taus[[jj]], function(x){length(which(dat[[jj]]$time == x))})
    Rs[[jj]] <- sapply(taus[[jj]], function(x){which(dat[[jj]]$time >= x)})
  }
  

  if(!is.null(nfolds)){
    cv.index <- list(NA)
    for(j in 1:J){
      fold1 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==1)))
      uncens <- which(dat[[j]]$status==1)
      cv.index[[j]] <- split(uncens, fold1)
      
      fold0 <- sample(rep(1:nfolds, length=sum(dat[[j]]$status==0)))
      cens <- which(dat[[j]]$status==0)
      temp0 <- split(cens, fold0)
      for(jj in 1:nfolds){
        cv.index[[j]][[jj]] <- c(cv.index[[j]][[jj]], temp0[[jj]])
      }
    }
  }
  
  
  # ------------------------------------------------------
  # Remove genes with no variability in at least one fold
  # ------------------------------------------------------
  if(!is.null(nfolds)){
    rm.genes <- NULL
    for(kk in 1:nfolds){
      for(j in 1:J){
        temp <- which(apply(dat[[j]]$X[cv.index[[j]][[kk]],], 2, sd) == 0)
        if(length(temp)> 0){
          rm.genes <- c(rm.genes, temp)
        }
      }
    }
  } else {
    rm.genes <- NULL
    for(j in 1:J){
      temp <- which(apply(dat[[j]]$X, 2, sd) == 0)
      if(length(temp)> 0){
        rm.genes <- c(rm.genes, temp)
      }
    }
  }
  
  if(length(rm.genes) > 0){
    for(j in 1:J){
      dat[[j]]$X <- dat[[j]]$X[,-rm.genes]
    }
  }
  p <- dim(dat[[1]]$X)[2]
  
  # ------------------------------------------------------
  # Standardize.x if necessary
  # ------------------------------------------------------
  if(standardize.x){
    # --- first check if need to 
    for(jj in 1:J){
      dat[[jj]]$x <- (dat[[jj]]$X - rep(1, dim(dat[[jj]]$X)[1])%*%t(colMeans(dat[[jj]]$X)))/(rep(1, dim(dat[[jj]]$X)[1])%*%t(apply(dat[[jj]]$X, 2, sd)))
    }
  } else {
    for(jj in 1:J){
      dat[[jj]]$x <- dat[[jj]]$X
    }
  }
  
  
  # ------------------------------------------------------
  # Get tuning candidate tuning parameter range 
  # ------------------------------------------------------
  grad.beta <- function(dat, beta, taus, Ds, Es, Rs, lambda){
    
    J <- length(dat)
    n <- sum(unlist(lapply(dat, function(x){dim(x$x)[1]})))
    gradout <- matrix(0, nrow=dim(beta)[1], ncol=dim(beta)[2])
    
    for(j in 1:J){
      
      t0 <- exp(tcrossprod(dat[[j]]$x, t(beta[,j])))
      t1 <- (dat[[j]]$x)*matrix(t0, nrow=length(t0), ncol=dim(dat[[j]]$x)[2])
      
      for(x in 1:length(taus[[j]])){
        temp0 <- sum(t0[Rs[[j]][[x]]])
        gradout[,j] <- gradout[,j] + (Ds[[j]][[x]]/temp0)*colsums(t1[Rs[[j]][[x]],,drop=FALSE]) - colsums(dat[[j]]$x[Es[[j]][[x]],,drop=FALSE])
      }
    }
    
    return(gradout/n)
  }
  
  ZetaList <- list(NA)
  beta.init <- matrix(0, nrow=p, ncol=J)
  for(j in 1:J){
    ZetaList[[j]] <- rep(0, dim(dat[[j]]$x)[1])
  }
  tempGrad <- grad.beta(dat, beta.init, taus, Ds, Es, Rs, lambda)

  lam.max <- max(apply(tempGrad, 1, function(x){sqrt(sum(x^2))}))
  lam.min <- delta2*lam.max
  lambda.vec <- 10^seq(log10(lam.max), log10(lam.min), length=nlambda2)
    
  
  # ------------------------------------------------------
  # Output matrices 
  # ------------------------------------------------------
  beta.out <- array(0, dim=c(1, length(lambda.vec), p, J))
  
  # ------------------------------------------------
  # Generate necessary preliminary quantities
  # ------------------------------------------------
  beta.temp <- beta.init

  for(kk in 1:length(lambda.vec)){
    
    temp <- beta.AccProxGrad.L2(dat = dat, beta.temp = beta.temp, taus = taus, Es = Es, 
                                Rs = Rs, Ds = Ds, lambda2 = lambda.vec[kk], 
                                max.iter = max.iter, tol = inner.tol, 
                                inner.quiet =   inner.quiet)
    beta.temp <- temp$beta
    beta.init <- beta.temp
    
    if(!quiet){
      cat("# ------------------------------------------------------", "\n")
      cat("Through tuning parameter", kk, " of ", length(lambda.vec),  "\n")
      cat(sum(rowSums(beta.temp!=0) > 0),":", sum(beta.temp^2), "\n")
    }
    
    beta.out[1,kk,,] <- beta.temp
    
  }

  
  if(!is.null(nfolds)){
    
    # -----------------------------------------
    # save metric for cross-validation 
    # -----------------------------------------
    PartLik <- array(NA, dim=c(length(lambda.vec), 1, nfolds, J))
    linPredTest <- list(NA)
    for(j in 1:J){
      linPredTest[[j]] <- array(0, dim = c(dim(dat[[j]]$X)[1], length(lambda.vec), 1))
    }
    linPredPartLik <- array(0, dim=c(length(lambda.vec), 1, J))
    
    for(kk in 1:nfolds){
    
      tausTrain <- list(NA)
      EsTrain <- list(NA)
      RsTrain <- list(NA)
      DsTrain <- list(NA)
      taus <- NULL
      Es <- NULL
      Rs <- NULL
      Ds <- NULL
      
      # ------------------------------------------------
      # Standardize and creating training/testing data
      # ------------------------------------------------
      for(jj in 1:J){
        
        tempTrain <- dat[[jj]]$X[-cv.index[[jj]][[kk]],]
        tempTest <- dat[[jj]]$X[cv.index[[jj]][[kk]],]
        if(!standardize.x){
          dat[[jj]]$x <- tempTrain
          dat[[jj]]$xTest <- tempTest
        } else {
          dat[[jj]]$x <- (tempTrain - rep(1, dim(tempTrain)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTrain)[1])%*%t(apply(tempTrain, 2, sd)))
          dat[[jj]]$xTest <- (tempTest - rep(1, dim(tempTest)[1])%*%t(colMeans(tempTrain)))/(rep(1, dim(tempTest)[1])%*%t(apply(tempTrain, 2, sd)))
        }
        dat[[jj]]$timeTrain <- dat[[jj]]$time[-cv.index[[jj]][[kk]]]
        dat[[jj]]$timeTest <- dat[[jj]]$time[cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTrain <- dat[[jj]]$status[-cv.index[[jj]][[kk]]]
        dat[[jj]]$statusTest <-  dat[[jj]]$status[cv.index[[jj]][[kk]]]
        tausTrain[[jj]] <- dat[[jj]]$timeTrain[which(dat[[jj]]$statusTrain==1)]
        EsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain == x)})
        DsTrain[[jj]] <- sapply(tausTrain[[jj]],function(x){length(which(dat[[jj]]$timeTrain == x))})
        RsTrain[[jj]] <- sapply(tausTrain[[jj]], function(x){which(dat[[jj]]$timeTrain >= x)})
        
      }
      
      # -----------------------------------------------
      # Get tuning candidate tuning parameter range 
      # -----------------------------------------------
      beta.temp <- matrix(0, nrow=p, ncol=J)

      for(ll in 1:length(lambda.vec)){
        
        temp0 <- beta.AccProxGrad.L2(dat = dat, beta.temp = beta.temp, taus = tausTrain, Es = EsTrain, 
                                     Rs = RsTrain, Ds = DsTrain, lambda2 = lambda.vec[ll], 
                                     max.iter = max.iter, tol = inner.tol, 
                                     inner.quiet = TRUE)
        beta.temp <- temp0$beta
        if(!inner.quiet){
          cat("Through ", ll, " of ", length(lambda.vec), " in fold ", kk, " of ", nfolds, "\n")
        }

        for(j in 1:J){
          dev1 <- coxnet.deviance(pred = dat[[j]]$x%*%temp0$beta[,j], y = Surv(dat[[j]]$timeTrain, dat[[j]]$statusTrain))
          dev2 <- coxnet.deviance(pred = dat[[j]]$X%*%temp0$beta[,j], y = Surv(dat[[j]]$time, dat[[j]]$status))
          PartLik[ll,1,  kk, j] <- dev2 - dev1
          linPredTest[[j]][cv.index[[j]][[kk]],ll,1] <- dat[[j]]$xTest%*%temp0$beta[,j]
        }
      }
      
      if(!quiet){
        cat("# ------------------------- ", "\n")
        cat("Through fold ", kk, " of ", nfolds, "\n")
      }
    }
    
    for(ll in 1:length(lambda.vec)){
      for(jjj in 1:J){
        linPredPartLik[ll,1,jjj] <- coxnet.deviance(y=Surv(dat[[jjj]]$time, dat[[jjj]]$status), pred=linPredTest[[jjj]][,ll,1])
      }
    }
    
    
  } else {
    PartLik <- NULL
    linPredPartLik <- NULL
    linPredTest <- NULL
  }
  
  
  
  
  
  return(list("PartLik" = PartLik, 
              "linPredTest" = linPredTest,
              "linPredPartLik" = linPredPartLik,
              "beta" = beta.out, 
              "rm.genes" = rm.genes,
              "X.mean" = lapply(dat, function(x){colMeans(x$X)}),
              "X.sd" = lapply(dat, function(x){apply(x$X, 2, sd)}),
              "lambda.vec" = 0, 
              "lambda2.vec" = lambda.vec, 
              "standardize" = standardize.x))
  
}





# temp1 <- RRCox.cv(dat, nlambda1 = 10, nlambda2 = 0,
#                   delta1 = .5, delta2 = NULL, lambda1.min = NULL,
#                   nfolds = 5, 
#                   standardize.x = TRUE,
#                   max.iter = 500,
#                   inner.quiet = FALSE, 
#                   quiet = FALSE,
#                   inner.tol = 1e-8)

# temp2 <- RRCox.cv(dat, nlambda1 = 0, nlambda2 = 10,
#                   delta1 = 0, delta2 = .5, lambda1.min = NULL,
#                   nfolds = 5, 
#                   standardize.x = TRUE,
#                   max.iter = 500,
#                   inner.quiet = FALSE, 
#                   quiet = FALSE,
#                   inner.tol = 1e-8)


# temp3 <- RRCox.cv(dat, nlambda1 = 4, nlambda2 = 10,
#                   delta1 = NULL, delta2 = .5, lambda1.min = 10^-4,
#                   nfolds = 5, 
#                   standardize.x = TRUE,
#                   max.iter = 500,
#                   inner.quiet = FALSE, 
#                   quiet = FALSE,
#                   inner.tol = 1e-8)


