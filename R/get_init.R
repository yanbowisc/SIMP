#' Warm-start estimator for MCMC algorithm of SIMP
#' 
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' @param X1D_ctr Centered design matrix of the discrete part of the predictors of interest.
#' @param X2_ctr Centered design matrix of the predictors of not main interest
#' @param Y Response data matrix.
#' @param dx Partial predictor envelope dimension. Must be an integer between 0 and p_C.
#' @param dy Partial response envelope dimension. Must be an integer between 0 and r.
#' @param method.idx Index for the warm start method, default to be 1.
#'
#' @export
get_init <- function(X1C, X1D_ctr, X2_ctr, Y, dx, dy, method.idx = 1){
  dim_Y <- dim(Y)
  n <- dim_Y[1]
  r <- dim_Y[2]
  pc <- ncol(X1C)
  if (is.null(X1D_ctr)){
    pd <- 0
  }else{
    pd <- ncol(X1D_ctr)
  }

  if (is.null(X2_ctr)){
    p2 <- 0
  }else{
    p2 <- ncol(X2_ctr)
  }

  X.order <- 1:pc
  Y.order <- 1:r
  A_exists <- (dx > 0) & (dx < pc)
  B_exists <- (dy > 0) & (dy < r)

  muX1C <- colMeans(X1C)
  X1Cc <- sweep(X1C, 2, muX1C, "-")
  X1c <- cbind(X1Cc, X1D_ctr)
  muY <- colMeans(Y)
  Yc <- sweep(Y, 2, muY, "-")

  if (pd > 0){
    XcXd.fit <- lm(X1Cc ~ X1D_ctr)
    gamma <- XcXd.fit$coefficients[-1,]
    row.names(gamma) <- NULL
    SigmaCD <- cov(XcXd.fit$residuals)
  }else{
    gamma <- 0
    SigmaCD <- cov(X1Cc)
  }

  if (p2 > 0){
    Yevlp.fit <- penv(X1c, X2_ctr, Yc, dy, asy = FALSE, init = NULL)
    beta2 <- t(Yevlp.fit$beta2)
    beta1C <- t(Yevlp.fit$beta1[,1:pc])
    if (pd > 0){
      beta1D <- t(Yevlp.fit$beta1[,(pc + 1):(pc + pd)])
    }else{
      beta1D <- 0
    }
    SigmaYX <- Yevlp.fit$Sigma
  }else{
    Yevlp.fit <- env(X1c, Yc, dy, asy = FALSE, init = NULL)
    beta1C <- t(Yevlp.fit$beta[,1:pc])
    if (pd > 0){
      beta1D <- t(Yevlp.fit$beta[,(pc + 1):(pc + pd)])
    }else{
      beta1D <- 0
    }
    beta2 <- 0
    SigmaYX <- Yevlp.fit$Sigma
  }

  if (B_exists){
    B_order <- find_A_from_gamma(Yevlp.fit$Gamma)
    B <- B_order$A
    if (B_order$switch.flag){
      Y.order <- B_order$Y.order
      Yc <- Yc[,Y.order]
      muY <- muY[Y.order]
      SigmaYX <- SigmaYX[Y.order, Y.order]
      if (p2 > 0){
        beta2 <- beta2[ ,Y.order]
      }
    }

    R_R0 <- find_gammas_from_A(B)
    R <- R_R0$gamma
    R0 <- R_R0$gamma0

    Phi <- crossprod(R, SigmaYX%*%R)
    Phi0 <- crossprod(R0, SigmaYX%*%R0)

    if (pd > 0){
      etaD <- t(beta1D%*%R)
    }else{
      etaD <- 0
    }

  }else if (dy == 0){
    B <- R <- Phi <- 0
    R0 <- diag(1, r)
    Phi0 <- SigmaYX
    etaD <- 0
  }else if (dy == r){
    B <- R0 <- Phi0 <- 0
    R <- diag(1, r)
    Phi <- SigmaYX
    if (pd > 0){
      etaD <- t(beta1D)
    }else{
      etaD <- 0
    }
  }

  if ((pd > 0)&(p2 > 0)){
    Xstar <- lm(X1Cc ~ X1D_ctr + X2_ctr)$residuals
  }else if ((pd == 0)&(p2 > 0)){
    Xstar <- lm(X1Cc ~ X2_ctr)$residuals
  }else if ((pd > 0)&(p2 == 0)){
    Xstar <- lm(X1Cc ~ X1D_ctr)$residuals
  }else if ((pd == 0)&(p2 == 0)){
    Xstar <- X1Cc
  }

  if (method.idx == 1){
    
    if ((pd > 0)&(p2 > 0)){
      Ystar <- lm(Yc ~ X1D_ctr + X2_ctr)$residuals
    }else if ((pd == 0)&(p2 > 0)){
      Ystar <- lm(Yc ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      Ystar <- lm(Yc ~ X1D_ctr)$residuals
    }else if ((pd == 0)&(p2 == 0)){
      Ystar <- Yc
    }
    Xenvlp.fit <- xenv(X = Xstar, Y = Ystar, dx, asy = FALSE, init = NULL)
    
  }else if(method.idx == 2){
    
    if ((pd > 0)&(p2 > 0)){
      YstarR <- lm(Yc%*%R ~ X1D_ctr + X2_ctr)$residuals
    }else if ((pd == 0)&(p2 > 0)){
      YstarR <- lm(Yc%*%R ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      YstarR <- lm(Yc%*%R ~ X1D_ctr)$residuals
    }else if ((pd == 0)&(p2 == 0)){
      YstarR <- Yc%*%R
    }
    Xenvlp.fit <- xenv(X = Xstar, Y = YstarR, dx, asy = FALSE, init = NULL)
    
  }else if(method.idx == 3){
    
    if ((pd > 0)&(p2 > 0)){
      Yevlp.fit1 <- penv(X1D_ctr, X2_ctr, Yc, dy, asy = FALSE, init = NULL)
      Ystar1 <- Yc - tcrossprod(X1D_ctr, Yevlp.fit1$beta1) - tcrossprod(X2_ctr, Yevlp.fit1$beta2)
    }else if ((pd == 0)&(p2 > 0)){
      Ystar1 <- lm(Yc ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      Yevlp.fit1 <- env(X1D_ctr, Yc, dy, asy = FALSE, init = NULL)
      Ystar1 <- Yc - tcrossprod(X1D_ctr, Yevlp.fit1$beta)
    }else if ((pd == 0)&(p2 == 0)){
      Ystar1 <- Yc
    }
    Xenvlp.fit <- xenv(X = Xstar, Y = Ystar1, dx, asy = FALSE, init = NULL)
    
  }else if(method.idx == 4){
    
    if ((pd > 0)&(p2 > 0)){
      Ystar <- lm(Yc ~ X1D_ctr + X2_ctr)$residuals
    }else if ((pd == 0)&(p2 > 0)){
      Ystar <- lm(Yc ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      Ystar <- lm(Yc ~ X1D_ctr)$residuals
    }else if ((pd == 0)&(p2 == 0)){
      Ystar <- Yc
    }
    Xenvlp.fit <- stenv(X = Xstar, Y = Ystar, dx, dy, asy = FALSE)
    
  }else if(method.idx == 5){
    
    if ((pd > 0)&(p2 > 0)){
      YstarR <- lm(Yc%*%R ~ X1D_ctr + X2_ctr)$residuals
    }else if ((pd == 0)&(p2 > 0)){
      YstarR <- lm(Yc%*%R ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      YstarR <- lm(Yc%*%R ~ X1D_ctr)$residuals
    }else if ((pd == 0)&(p2 == 0)){
      YstarR <- Yc%*%R
    }
    Xenvlp.fit <- stenv(X = Xstar, Y = YstarR, dx, dy, asy = FALSE)
    
  }else if(method.idx == 6){
    
    if ((pd > 0)&(p2 > 0)){
      Yevlp.fit1 <- penv(X1D_ctr, X2_ctr, Yc, dy, asy = FALSE, init = NULL)
      Ystar1 <- Yc - tcrossprod(X1D_ctr, Yevlp.fit1$beta1) - tcrossprod(X2_ctr, Yevlp.fit1$beta2)
    }else if ((pd == 0)&(p2 > 0)){
      Ystar1 <- lm(Yc ~ X2_ctr)$residuals
    }else if ((pd > 0)&(p2 == 0)){
      Yevlp.fit1 <- env(X1D_ctr, Yc, dy, asy = FALSE, init = NULL)
      Ystar1 <- Yc - tcrossprod(X1D_ctr, Yevlp.fit1$beta)
    }else if ((pd == 0)&(p2 == 0)){
      Ystar1 <- Yc
    }
    Xenvlp.fit <- stenv(X = Xstar, Y = Ystar1, dx, dy, asy = FALSE)
    
  }


  if (A_exists){
    if (method.idx %in% 1:3){
      A_order <- find_A_from_gamma(Xenvlp.fit$Gamma)
    }else if(method.idx %in% 4:6){
      A_order <- find_A_from_gamma(Xenvlp.fit$Phi)
    }

    A <- A_order$A

    if (A_order$switch.flag){
      X.order <- A_order$Y.order
      X1Cc <- X1Cc[ ,X.order]
      muX1C <- muX1C[X.order]
      gamma <- gamma[ ,X.order]
      SigmaCD <- SigmaCD[X.order, X.order]
    }

    L_L0 <- find_gammas_from_A(A)
    L <- L_L0$gamma
    L0 <- L_L0$gamma0

    Omega <- crossprod(L, SigmaCD%*%L)
    Omega0 <- crossprod(L0, SigmaCD%*%L0)

    if (dy > 0){
      etaC <- crossprod(beta1C%*%R, L)
    }else{
      etaC <- 0
    }

  }else if (dx == 0){
    A <- L <- Omega <- 0
    L0 <- diag(1, pc)
    Omega0 <- SigmaCD
    etaC <- 0
  }else if (dx == pc){
    A <- L0 <- Omega0 <- 0
    L <- diag(1, pc)
    Omega <- SigmaCD

    if (dy > 0){
      etaC <- crossprod(beta1C%*%R, L)
    }else{
      etaC <- 0
    }
  }

  list(muX1C = muX1C,
       muY = muY,
       beta1C = beta1C,
       beta1D = beta1D,
       beta2 = beta2,
       gamma = gamma,
       Omega = Omega,
       Omega0 = Omega0,
       Phi = Phi,
       Phi0 = Phi0,
       A = A,
       L = L,
       L0 = L0,
       B = B,
       R = R,
       R0 = R0,
       etaC = etaC,
       etaD = etaD,
       SigmaCD = SigmaCD,
       SigmaYX = SigmaYX,
       X1Cc = X1Cc,
       Yc = Yc,
       X.order = X.order,
       Y.order = Y.order)

}
