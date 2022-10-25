#' BIC calculation for the output of \code{SIMP()}
#' 
#' @param SIMP.fit The output of \code{SIMP()}.
#' @param r Dimension of response Y.
#' @param pc Dimension of X1C, the continuous part of the predictors of interest.
#' @param pd Dimension of X1D, the discrete part of the predictors of interest.
#' @param p2 Dimension of X2, predictors of not main interest.
#' @param dx Partial predictor envelope dimension of SIMP.
#' @param dy Partial response envelope dimension of SIMP.
#' @param samp.size Sample size of the datasets.
#' @param combined Logical. Indicate whether posterior samples (after burn-in) from all chains 
#' should be combined to calculate BIC. If there is only one chain, then there is no difference 
#' by indicating this parameter to be TRUE or FALSE.
#' @param Chain.no If combined = FALSE, we need to indicate which chain should be used for the 

#' calculation of BIC-MCMC. Chain.no must be an integer between 1 and the number of chains.
#'
#' @export
BIC <- function(SIMP.fit, r, pc, pd, p2, dx, dy, samp.size, combined = T, Chain.no = 1){
  
  if (combined){
    llik.combined <- do.call(c, SIMP.fit$llik)
    names(llik.combined) <- NULL
  }else{
    llik.combined <- SIMP.fit$llik[[Chain.no]]
  }
  k <- dy * (dx + pd) + r * (r + 2 * p2 + 3)/2 + pc * (pc + 3)/2
  bic <- -2 * max(llik.combined) + k * log(samp.size)
  bic
}

#' AIC-MCMC calculation for the output of \code{SIMP()}
#' 
#' @param SIMP.fit The output of \code{SIMP()}.
#' @param r Dimension of response Y.
#' @param pc Dimension of X1C, the continuous part of the predictors of interest.
#' @param pd Dimension of X1D, the discrete part of the predictors of interest.
#' @param p2 Dimension of X2, predictors of not main interest.
#' @param dx Partial predictor envelope dimension of SIMP.
#' @param dy Partial response envelope dimension of SIMP.
#' @param combined Logical. Indicate whether posterior samples (after burn-in) from all chains 
#' should be combined to calculate AIC. If there is only one chain, then there is no difference 
#' by indicating this parameter to be TRUE or FALSE.
#' @param Chain.no If combined = FALSE, we need to indicate which chain should be used for the 
#' calculation of AIC. Chain.no must be an integer between 1 and the number of chains.
#'
#' @export
AIC <- function(SIMP.fit, r, pc, pd, p2, dx, dy, combined = T, Chain.no = 1){
  if (combined){
    llik.combined <- do.call(c, SIMP.fit$llik)
    names(llik.combined) <- NULL
  }else{
    llik.combined <- SIMP.fit$llik[[Chain.no]]
  }
  k <- dy * (dx + pd) + r * (r + 2 * p2 + 3)/2 + pc * (pc + 3)/2
  aic <- -2 * max(llik.combined) + 2 * k
  aic
}

#' DIC calculation for the output of \code{SIMP()}
#' 
#' @param SIMP.fit The output of \code{SIMP()}.
#'
#' @param Y Response matrix. Must have the same number of rows as \code{X}.
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X1D Design matrix of the discrete part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X2 Design matrix of the nuisance predictors. Must have the same number of rows as \code{Y}.
#' @param dx Partial predictor envelope dimension. Must be an integer between 0 and \code{ncol(X1C)}.
#' @param dy Partial response envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' @param combined Logical. Indicate whether posterior samples (after burn-in) from all chains 
#' should be combined to calculate DIC. If there is only one chain, then there is no difference 
#' by indicating this parameter to be TRUE or FALSE.
#' @param Chain.no If combined = FALSE, we need to indicate which chain should be used for the 
#' calculation of DIC. Chain.no must be an integer between 1 and the number of chains.
#'
#' @export
DIC <- function(SIMP.fit, Y, X1C, X1D, X2, dx, dy, combined = T, Chain.no = 1){
  n <- nrow(Y)
  pc <- ncol(X1C)
  r <- ncol(Y)
  if(is.null(X1D)){
    pd <- 0
  }else{
    pd <- ncol(X1D)
    X1D_ctr <- scale(X1D, center = TRUE, scale = FALSE)
  }
  
  if(is.null(X2)){
    p2 <- 0
  }else{
    p2 <- ncol(X2)
    X2_ctr <- scale(X2, center = TRUE, scale = FALSE)
  }
  
  if (combined){
    Omega.combined <- do.call(c, SIMP.fit$Omega)
    names(Omega.combined) <- NULL
  }else{
    Omega.combined <- SIMP.fit$Omega[[Chain.no]]
  }
  
  if (combined){
    Omega0.combined <- do.call(c, SIMP.fit$Omega0)
    names(Omega0.combined) <- NULL
  }else{
    Omega0.combined <- SIMP.fit$Omega0[[Chain.no]]
  }
  
  if (combined){
    Phi.combined <- do.call(c, SIMP.fit$Phi)
    names(Phi.combined) <- NULL
  }else{
    Phi.combined <- SIMP.fit$Phi[[Chain.no]]
  }
  
  if (combined){
    Phi0.combined <- do.call(c, SIMP.fit$Phi0)
    names(Phi0.combined) <- NULL
  }else{
    Phi0.combined <- SIMP.fit$Phi0[[Chain.no]]
  }
  
  Omega_est <- elmwise_mean_in_list(Omega.combined)
  Omega0_est <- elmwise_mean_in_list(Omega0.combined)
  Phi_est <- elmwise_mean_in_list(Phi.combined)
  Phi0_est <- elmwise_mean_in_list(Phi0.combined)
  
  log_det_Omega_est <- ifelse(dx > 0, log(det(Omega_est)), 0)
  log_det_Omega0_est <- ifelse(dx < pc, log(det(Omega0_est)), 0)
  log_det_Phi_est <- ifelse(dy > 0, log(det(Phi_est)), 0)
  log_det_Phi0_est <- ifelse(dy < r, log(det(Phi0_est)), 0)
  
  if (combined){
    muX1C.combined <- do.call(c, SIMP.fit$muX1C)
    names(muX1C.combined) <- NULL
  }else{
    muX1C.combined <- SIMP.fit$muX1C[[Chain.no]]
  }
  
  X1C_ctr <- as.matrix(sweep(X1C, 2, elmwise_mean_in_list(muX1C.combined), "-"))
  
  if (combined){
    gamma.combined <- do.call(c, SIMP.fit$gamma)
    names(gamma.combined) <- NULL
  }else{
    gamma.combined <- SIMP.fit$gamma[[Chain.no]]
  }
  gamma_est <- elmwise_mean_in_list(gamma.combined)
  
  if (pd > 0){
    X1C_ctr_shifted <- X1C_ctr - X1D_ctr%*%gamma_est
  }else{
    X1C_ctr_shifted <- X1C_ctr
  }
  
  if (combined){
    SigmaCD.combined <- do.call(c, SIMP.fit$SigmaCD)
    names(SigmaCD.combined) <- NULL
  }else{
    SigmaCD.combined <- SIMP.fit$SigmaCD[[Chain.no]]
  }
  
  if (combined){
    SigmaYX.combined <- do.call(c, SIMP.fit$SigmaYX)
    names(SigmaYX.combined) <- NULL
  }else{
    SigmaYX.combined <- SIMP.fit$SigmaYX[[Chain.no]]
  }
  
  SigmaCD.inv_est <- solve_chol(elmwise_mean_in_list(SigmaCD.combined))
  SigmaYX.inv_est <- solve_chol(elmwise_mean_in_list(SigmaYX.combined))
  
  if (combined){
    muY.combined <- do.call(c, SIMP.fit$muY)
    names(muY.combined) <- NULL
  }else{
    muY.combined <- SIMP.fit$muY[[Chain.no]]
  }
  Yctr <- as.matrix(sweep(Y, 2, elmwise_mean_in_list(muY.combined), "-"))
  
  
  if (combined){
    beta1C.combined <- do.call(c, SIMP.fit$beta1C)
    names(beta1C.combined) <- NULL
  }else{
    beta1C.combined <- SIMP.fit$beta1C[[Chain.no]]
  }
  
  if (combined){
    beta1D.combined <- do.call(c, SIMP.fit$beta1D)
    names(beta1D.combined) <- NULL
  }else{
    beta1D.combined <- SIMP.fit$beta1D[[Chain.no]]
  }
  
  if (combined){
    beta2.combined <- do.call(c, SIMP.fit$beta2)
    names(beta2.combined) <- NULL
  }else{
    beta2.combined <- SIMP.fit$beta2[[Chain.no]]
  }

  beta1C_est <- elmwise_mean_in_list(beta1C.combined)
  beta1D_est <- elmwise_mean_in_list(beta1D.combined)
  beta2_est <- elmwise_mean_in_list(beta2.combined)
  
  if (pd * p2 > 0){
    resi <- Yctr - X1C_ctr%*%beta1C_est - X1D_ctr%*%beta1D_est - X2_ctr%*%beta2_est
  }else if ((pd == 0)&(p2 > 0)){
    resi <- Yctr - X1C_ctr%*%beta1C_est - X2_ctr%*%beta2_est
  }else if ((pd > 0)&(p2 == 0)){
    resi <- Yctr - X1C_ctr%*%beta1C_est - X1D_ctr%*%beta1D_est
  }else{
    resi <- Yctr - X1C_ctr%*%beta1C_est
  }
  
  llik_bayes <- -0.5*(n*(log_det_Omega_est + log_det_Omega0_est) +
                sum(X1C_ctr_shifted %*% SigmaCD.inv_est * X1C_ctr_shifted) +
                n*(log_det_Phi_est + log_det_Phi0_est) +
                sum(resi %*% SigmaYX.inv_est * resi))
  
  if (combined){
    llik.combined <- do.call(c, SIMP.fit$llik)
    names(llik.combined) <- NULL
  }else{
    llik.combined <- SIMP.fit$llik[[Chain.no]]
  }
  
  p_DIC <- 2*(llik_bayes - mean(llik.combined))
  dic <- -2 * llik_bayes + 2 * p_DIC
  dic
}

#' WAIC calculation for the output of \code{SIMP()}
#' 
#' @param SIMP.fit The output of \code{SIMP()}.
#'
#' @param Y Response matrix. Must have the same number of rows as \code{X}.
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X1D Design matrix of the discrete part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X2 Design matrix of the nuisance predictors. Must have the same number of rows as \code{Y}.
#' @param dx Partial predictor envelope dimension. Must be an integer between 0 and \code{ncol(X1C)}.
#' @param dy Partial response envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' @param combined Logical. Indicate whether posterior samples (after burn-in) from all chains 
#' should be combined to calculate WAIC. If there is only one chain, then there is no difference 
#' by indicating this parameter to be TRUE or FALSE.
#' @param Chain.no If combined = FALSE, we need to indicate which chain should be used for the 
#' calculation of WAIC. Chain.no must be an integer between 1 and the number of chains.
#'
#' @export
WAIC <- function(SIMP.fit, Y, X1C, X1D, X2, dx, dy, combined = T, Chain.no = 1){
  Y <- as.matrix(Y)
  X1C <- as.matrix(X1C)
  samp.size <- nrow(Y)
  if (combined){
    muY.combined <- do.call(c, SIMP.fit$muY)
    names(muY.combined) <- NULL
  }else{
    muY.combined <- SIMP.fit$muY[[Chain.no]]
  }
  S <- length(muY.combined)
  pc <- ncol(X1C)
  r <- ncol(Y)
  if(is.null(X1D)){
    pd <- 0
  }else{
    pd <- ncol(X1D)
    X1D_ctr <- scale(X1D, center = TRUE, scale = FALSE)
  }
  
  if(is.null(X2)){
    p2 <- 0
  }else{
    p2 <- ncol(X2)
    X2_ctr <- scale(X2, center = TRUE, scale = FALSE)
  }
  
  if (dx > 0){
    if (combined){
      Omega.combined <- do.call(c, SIMP.fit$Omega)
      names(Omega.combined) <- NULL
    }else{
      Omega.combined <- SIMP.fit$Omega[[Chain.no]]
    }
    log_det_Omega <- log(unlist(lapply(Omega.combined, det)))
  }else{
    log_det_Omega <- rep(0,  S)
  }
  
  if (dx < pc){
    if (combined){
      Omega0.combined <- do.call(c, SIMP.fit$Omega0)
      names(Omega0.combined) <- NULL
    }else{
      Omega0.combined <- SIMP.fit$Omega0[[Chain.no]]
    }
    log_det_Omega0 <- log(unlist(lapply(Omega0.combined, det)))
  }else{
    log_det_Omega0 <- rep(0,  S)
  }
  
  if (dy > 0){
    if (combined){
      Phi.combined <- do.call(c, SIMP.fit$Phi)
      names(Phi.combined) <- NULL
    }else{
      Phi.combined <- SIMP.fit$Phi[[Chain.no]]
    }
    log_det_Phi <- log(unlist(lapply(Phi.combined, det)))
  }else{
    log_det_Phi <- rep(0,  S)
  }
  
  if (dy < r){
    if (combined){
      Phi0.combined <- do.call(c, SIMP.fit$Phi0)
      names(Phi0.combined) <- NULL
    }else{
      Phi0.combined <- SIMP.fit$Phi0[[Chain.no]]
    }
    log_det_Phi0 <- log(unlist(lapply(Phi0.combined, det)))
  }else{
    log_det_Phi0 <- rep(0,  S)
  }
  
  if (combined){
    SigmaCD.combined <- do.call(c, SIMP.fit$SigmaCD)
    names(SigmaCD.combined) <- NULL
  }else{
    SigmaCD.combined <- SIMP.fit$SigmaCD[[Chain.no]]
  }
  
  if (combined){
    SigmaYX.combined <- do.call(c, SIMP.fit$SigmaYX)
    names(SigmaYX.combined) <- NULL
  }else{
    SigmaYX.combined <- SIMP.fit$SigmaYX[[Chain.no]]
  }
  
  SigmaCD.inv <- lapply(SigmaCD.combined, solve_chol)
  SigmaYX.inv <- lapply(SigmaYX.combined, solve_chol)
  
  if (combined){
    gamma.combined <- do.call(c, SIMP.fit$gamma)
    names(gamma.combined) <- NULL
  }else{
    gamma.combined <- SIMP.fit$gamma[[Chain.no]]
  }
  
  if (combined){
    beta1C.combined <- do.call(c, SIMP.fit$beta1C)
    names(beta1C.combined) <- NULL
  }else{
    beta1C.combined <- SIMP.fit$beta1C[[Chain.no]]
  }
  
  if (combined){
    beta1D.combined <- do.call(c, SIMP.fit$beta1D)
    names(beta1D.combined) <- NULL
  }else{
    beta1D.combined <- SIMP.fit$beta1D[[Chain.no]]
  }
  
  if (combined){
    beta2.combined <- do.call(c, SIMP.fit$beta2)
    names(beta2.combined) <- NULL
  }else{
    beta2.combined <- SIMP.fit$beta2[[Chain.no]]
  }
  
  if (combined){
    muX1C.combined <- do.call(c, SIMP.fit$muX1C)
    names(muX1C.combined) <- NULL
  }else{
    muX1C.combined <- SIMP.fit$muX1C[[Chain.no]]
  }
  
  lppd <- 0
  p_WAIC <- 0
  for (i in 1:samp.size){
    llik_i <- numeric(S)
    for (s in 1:S){
      
      X1C_ctr_i_s <- X1C[i, ] - muX1C.combined[[s]]
      
      if (pd > 0){
        X1C_ctr_shifted_i_s <- X1C_ctr_i_s - X1D_ctr[i, ]%*%gamma.combined[[s]]
      }else{
        X1C_ctr_shifted_i_s <- X1C_ctr_i_s
      }
      
      if (pd * p2 > 0){
        resi_i_s <- Y[i, ] - muY.combined[[s]]  - X1C_ctr_i_s%*%beta1C.combined[[s]] - X1D_ctr[i,]%*%beta1D.combined[[s]] - X2_ctr[i,]%*%beta2.combined[[s]]
      }else if ((pd == 0)&(p2 > 0)){
        resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]] - X2_ctr[i,]%*%beta2.combined[[s]]
      }else if ((pd > 0)&(p2 == 0)){
        resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]] - X1D_ctr[i,]%*%beta1D.combined[[s]]
      }else{
        resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]]
      }
      
      llik_i[s] <- - 0.5*(log_det_Omega[s] + log_det_Omega0[s] + log_det_Phi[s] + log_det_Phi0[s] + sum(X1C_ctr_shifted_i_s %*% SigmaCD.inv[[s]] * X1C_ctr_shifted_i_s) + sum(resi_i_s %*% SigmaYX.inv[[s]] * resi_i_s) )
    }
    lppd <- lppd + log(mean(exp(llik_i)))
    p_WAIC <- p_WAIC + var(llik_i)
  }
  
  WAIC <- -2 * lppd + 2 * p_WAIC
  WAIC
}

#' Bayesian CV 
#' @param SIMP.fit.full The output of \code{SIMP()} fitted on the full dataset. This input is used as
#' correction, and the posterior samples from all chains will be combined.
#'
#' @param Y Response matrix. Must have the same number of rows as \code{X}.
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X1D Design matrix of the discrete part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X2 Design matrix of the nuisance predictors. Must have the same number of rows as \code{Y}.
#' @param dx Partial predictor envelope dimension. Must be an integer between 0 and \code{ncol(X1C)}.
#' @param dy Partial response envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' @param K_CV Fold for CV
#' @param n.iter Number of Markov chain iterations to run in each chains, for the fitting in CV.
#'  *Includes burn-in*.
#' @param correction Whether correction should be applied.
#'
#' @export
CV <- function(SIMP.fit.full, Y, X1C, X1D, X2, dx, dy, K_CV, n.iter, correction = FALSE){
  Y <- as.matrix(Y)
  X1C <- as.matrix(X1C)
  X1D <- as.matrix(X1D)
  X2 <- as.matrix(X2)
  samp.size <- nrow(Y)
  pc <- ncol(X1C)
  r <- ncol(Y)
  lpd <- numeric(K_CV)
  lppd_i <- numeric(K_CV)
  for (k in 1:K_CV){
    
    X1C.sub <- X1C[-((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
    
    if (!is.null(X1D)){
      pd <- ncol(X1D)
      X1D.sub <- as.matrix(X1D[-((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ])
    }else{
      pd <- 0
      X1D.sub <- NULL
    }
    
    if (!is.null(X2)){
      p2 <- ncol(X2)
      X2.sub <- X2[-((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
    }else{
      p2 <- 0
      X2.sub <- NULL
    }
    
    Y.sub <- Y[-((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
    SIMP.fit.sub <- SIMP(X1C.sub, 
                         X1D.sub, 
                         X2.sub, 
                         Y.sub, 
                         dx = dx, 
                         dy = dy,
                         n.iter = n.iter,
                         n.chains = 1,
                         tau = 0.1,
                         init_method = "envlps",
                         Metro_method = "RW",
                         HMC_steps = NA,
                         autotune = TRUE,
                         tune.accpt.prop.lower = 0.3,
                         tune.accpt.prop.upper = 0.4,
                         tune.incr = 0.05,
                         burnin.prop = 0.5,
                         tune.burnin.prop = 0.5,
                         tune_nterm = 50,
                         show_progress = TRUE,
                         chains_parallel = FALSE)
    
    S <- length(SIMP.fit.sub$muY[[1]])
    
    if (dx > 0){
      log_det_Omega <- log(unlist(lapply(SIMP.fit.sub$Omega[[1]], det)))
    }else{
      log_det_Omega <- rep(0,  S)
    }
    
    if (dx < pc){
      log_det_Omega0 <- log(unlist(lapply(SIMP.fit.sub$Omega0[[1]], det)))
    }else{
      log_det_Omega0 <- rep(0,  S)
    }
    
    if (dy > 0){
      log_det_Phi <- log(unlist(lapply(SIMP.fit.sub$Phi[[1]], det)))
    }else{
      log_det_Phi <- rep(0,  S)
    }
    
    if (dy < r){
      log_det_Phi0 <- log(unlist(lapply(SIMP.fit.sub$Phi0[[1]], det)))
    }else{
      log_det_Phi0 <- rep(0,  S)
    }

    SigmaCD.inv <- lapply(SIMP.fit.sub$SigmaCD[[1]], solve_chol)
    SigmaYX.inv <- lapply(SIMP.fit.sub$SigmaYX[[1]], solve_chol)
    
    holdout_samp.size <- length((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV))
    
    X1C.sub.out <- X1C[((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
    
    if (!is.null(X1D)){
      X1D.sub.out <- X1D[((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
      X1D.sub.out_ctr <- scale(X1D.sub.out, center = TRUE, scale = FALSE)
    }else{
      X1D.sub.out <- NULL
    }
    
    if (!is.null(X2)){
      X2.sub.out <- X2[((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
      X2.sub.out_ctr <- scale(X2.sub.out, center = TRUE, scale = FALSE)
    }else{
      X2.sub.out <- NULL
    }
    
    Y.sub.out <- Y[((floor((k-1)*samp.size/K_CV) + 1):floor(k*samp.size/K_CV)), ]
    
    llik <- numeric(S)
      for (s in 1:S){
        
        X1C.sub.out_ctr_s <- X1C.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muX1C[[1]][[s]])
        
        if (pd > 0){
          X1C.sub.out_ctr_shifted_s <- X1C.sub.out_ctr_s - X1D.sub.out_ctr%*%SIMP.fit.sub$gamma[[1]][[s]]
        }else{
          X1C.sub.out_ctr_shifted_s <- X1C.sub.out_ctr_s
        }
        
        if (pd * p2 > 0){
          resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X1D.sub.out_ctr%*%SIMP.fit.sub$beta1D[[1]][[s]] - X2.sub.out_ctr%*%SIMP.fit.sub$beta2[[1]][[s]]
        }else if ((pd == 0)&(p2 > 0)){
          resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X2.sub.out_ctr%*%SIMP.fit.sub$beta2[[1]][[s]]
        }else if ((pd > 0)&(p2 == 0)){
          resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X1D.sub.out_ctr%*%SIMP.fit.sub$beta1D[[1]][[s]]
        }else{
          resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]]
        }
        
        llik[s] <- -0.5*(holdout_samp.size*(log_det_Omega[s] + log_det_Omega0[s]) +
                        sum(X1C.sub.out_ctr_shifted_s %*% SigmaCD.inv[[s]] * X1C.sub.out_ctr_shifted_s) +
                        holdout_samp.size*(log_det_Phi[s] + log_det_Phi0[s]) +
                        sum(resi.sub.out_s %*% SigmaYX.inv[[s]] * resi.sub.out_s))
      }
      lpd[k] <- log(mean(exp(llik - max(llik)))) + max(llik)
      
      if (correction == TRUE){
        
        lppd_i_each <- numeric(K_CV)
        
        for (l in 1:K_CV){
          holdout_samp.size <- length((floor((l-1)*samp.size/K_CV) + 1):floor(l*samp.size/K_CV))
          
          X1C.sub.out <- X1C[((floor((l-1)*samp.size/K_CV) + 1):floor(l*samp.size/K_CV)), ]
          
          if (!is.null(X1D)){
            X1D.sub.out <- X1D[((floor((l-1)*samp.size/K_CV) + 1):floor(l*samp.size/K_CV)), ]
            X1D.sub.out_ctr <- scale(X1D.sub.out, center = TRUE, scale = FALSE)
          }else{
            X1D.sub.out <- NULL
          }
          
          if (!is.null(X2)){
            X2.sub.out <- X2[((floor((l-1)*samp.size/K_CV) + 1):floor(l*samp.size/K_CV)), ]
            X2.sub.out_ctr <- scale(X2.sub.out, center = TRUE, scale = FALSE)
          }else{
            X2.sub.out <- NULL
          }
          
          Y.sub.out <- Y[((floor((l-1)*samp.size/K_CV) + 1):floor(l*samp.size/K_CV)), ]
          
          llik <- numeric(S)
          for (s in 1:S){
            
            X1C.sub.out_ctr_s <- X1C.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muX1C[[1]][[s]])
            
            if (pd > 0){
              X1C.sub.out_ctr_shifted_s <- X1C.sub.out_ctr_s - X1D.sub.out_ctr%*%SIMP.fit.sub$gamma[[1]][[s]]
            }else{
              X1C.sub.out_ctr_shifted_s <- X1C.sub.out_ctr_s
            }
            
            if (pd * p2 > 0){
              resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X1D.sub.out_ctr%*%SIMP.fit.sub$beta1D[[1]][[s]] - X2.sub.out_ctr%*%SIMP.fit.sub$beta2[[1]][[s]]
            }else if ((pd == 0)&(p2 > 0)){
              resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X2.sub.out_ctr%*%SIMP.fit.sub$beta2[[1]][[s]]
            }else if ((pd > 0)&(p2 == 0)){
              resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]] - X1D.sub.out_ctr%*%SIMP.fit.sub$beta1D[[1]][[s]]
            }else{
              resi.sub.out_s <- Y.sub.out - tcrossprod(rep(1, holdout_samp.size), SIMP.fit.sub$muY[[1]][[s]]) - X1C.sub.out_ctr_s%*%SIMP.fit.sub$beta1C[[1]][[s]]
            }
            
            llik[s] <- -0.5*(holdout_samp.size*(log_det_Omega[s] + log_det_Omega0[s]) +
                               sum(X1C.sub.out_ctr_shifted_s %*% SigmaCD.inv[[s]] * X1C.sub.out_ctr_shifted_s) +
                               holdout_samp.size*(log_det_Phi[s] + log_det_Phi0[s]) +
                               sum(resi.sub.out_s %*% SigmaYX.inv[[s]] * resi.sub.out_s))
          }
          lppd_i_each[l] <- log(mean(exp(llik - max(llik)))) + max(llik)
          
        }
        lppd_i[k] <- mean(lppd_i_each) 

    
        
      }
  }
  
  cv <- - sum(lpd)
  
  if (correction == TRUE){
    muY.combined <- do.call(c, SIMP.fit.full$muY)
    S <- length(muY.combined)
    if(is.null(X1D)){
      pd <- 0
    }else{
      pd <- ncol(X1D)
      X1D_ctr <- scale(X1D, center = TRUE, scale = FALSE)
    }
    
    if(is.null(X2)){
      p2 <- 0
    }else{
      p2 <- ncol(X2)
      X2_ctr <- scale(X2, center = TRUE, scale = FALSE)
    }
    
    if (dx > 0){
      Omega.combined <- do.call(c, SIMP.fit.full$Omega)
      names(Omega.combined) <- NULL
      log_det_Omega <- log(unlist(lapply(Omega.combined, det)))
    }else{
      log_det_Omega <- rep(0,  S)
    }
    
    if (dx < pc){
      Omega0.combined <- do.call(c, SIMP.fit.full$Omega0)
      names(Omega0.combined) <- NULL
      log_det_Omega0 <- log(unlist(lapply(Omega0.combined, det)))
    }else{
      log_det_Omega0 <- rep(0,  S)
    }
    
    if (dy > 0){
      Phi.combined <- do.call(c, SIMP.fit.full$Phi)
      names(Phi.combined) <- NULL
      log_det_Phi <- log(unlist(lapply(Phi.combined, det)))
    }else{
      log_det_Phi <- rep(0,  S)
    }
    
    if (dy < r){
      Phi0.combined <- do.call(c, SIMP.fit.full$Phi0)
      names(Phi0.combined) <- NULL
      log_det_Phi0 <- log(unlist(lapply(Phi0.combined, det)))
    }else{
      log_det_Phi0 <- rep(0,  S)
    }
    
    SigmaCD.combined <- do.call(c, SIMP.fit.full$SigmaCD)
    names(SigmaCD.combined) <- NULL
    SigmaCD.inv <- lapply(SigmaCD.combined, solve_chol)
    
    SigmaYX.combined <- do.call(c, SIMP.fit.full$SigmaYX)
    names(SigmaYX.combined) <- NULL
    SigmaYX.inv <- lapply(SigmaYX.combined, solve_chol)
    
    muX1C.combined <- do.call(c, SIMP.fit.full$muX1C)
    names(muX1C.combined) <- NULL
    
    gamma.combined <- do.call(c, SIMP.fit.full$gamma)
    names(gamma.combined) <- NULL
    
    beta1C.combined <- do.call(c, SIMP.fit.full$beta1C)
    names(beta1C.combined) <- NULL
    
    beta1D.combined <- do.call(c, SIMP.fit.full$beta1D)
    names(beta1D.combined) <- NULL
    
    beta2.combined <- do.call(c, SIMP.fit.full$beta2)
    names(beta2.combined) <- NULL
    
    lppd <- 0
    for (i in 1:samp.size){
      llik_i <- numeric(S)
      
      
      for (s in 1:S){
        
        X1C_ctr_i_s <- X1C[i, ] - muX1C.combined[[s]]
        
        if (pd > 0){
          X1C_ctr_shifted_i_s <- X1C_ctr_i_s - X1D_ctr[i, ]%*%gamma.combined[[s]]
        }else{
          X1C_ctr_shifted_i_s <- X1C_ctr_i_s
        }
        
        if (pd * p2 > 0){
          resi_i_s <- Y[i, ] - muY.combined[[s]]  - X1C_ctr_i_s%*%beta1C.combined[[s]] - X1D_ctr[i,]%*%beta1D.combined[[s]] - X2_ctr[i,]%*%beta2.combined[[s]]
        }else if ((pd == 0)&(p2 > 0)){
          resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]] - X2_ctr[i,]%*%beta2.combined[[s]]
        }else if ((pd > 0)&(p2 == 0)){
          resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]] - X1D_ctr[i,]%*%beta1D.combined[[s]]
        }else{
          resi_i_s <- Y[i, ] - muY.combined[[s]] - X1C_ctr_i_s%*%beta1C.combined[[s]]
        }
        
        llik_i[s] <- - 0.5*(log_det_Omega[s] + log_det_Omega0[s] + log_det_Phi[s] + log_det_Phi0[s] + sum(X1C_ctr_shifted_i_s %*% SigmaCD.inv[[s]] * X1C_ctr_shifted_i_s) + sum(resi_i_s %*% SigmaYX.inv[[s]] * resi_i_s) )
      }
      lppd <- lppd + log(mean(exp(llik_i)))
    }
    
    
    b <- lppd - sum(lppd_i)
    cv <- cv - b
    
  }
  
  cv
}



