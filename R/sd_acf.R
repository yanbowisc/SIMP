#' Function to calculate the standard deviations and the estimates for the auto-correlation function for 
#' posterior samples (PSD and PACF).
#' @param SIMP.fit Returned outputs from SIMP().
#' @param thinning_length The length for thinning before calculating PSD and PACF.
#' @param acf.lagmax The maximal lag length for calculating PACF.
#'
#' @export
sd_acf <- function(SIMP.fit, thinning_length=1, acf.lagmax = 30){
  
  pos.samp.size <- length(SIMP.fit$muX1C)
  idx.after.thinned <- seq(1, pos.samp.size, by = thinning_length)
  SIMP.fit$muX1C <- SIMP.fit$muX1C[idx.after.thinned]
  SIMP.fit$muY <- SIMP.fit$muY[idx.after.thinned]
  SIMP.fit$beta1C <- SIMP.fit$beta1C[idx.after.thinned]
  SIMP.fit$beta1D <- SIMP.fit$beta1D[idx.after.thinned]
  SIMP.fit$beta2 <- SIMP.fit$beta2[idx.after.thinned]
  SIMP.fit$SigmaCD <- SIMP.fit$SigmaCD[idx.after.thinned]
  SIMP.fit$SigmaYX <- SIMP.fit$SigmaYX[idx.after.thinned]
  
  
  muX1C.flat <- matrix(unlist(SIMP.fit$muX1C), ncol = length(SIMP.fit$muX1C[[1]]), byrow = T)
  muX1C.sd <- apply(muX1C.flat, 2, sd)
  muX1C.acf <- apply(muX1C.flat, 2, function(x) acf(x, lag.max = acf.lagmax))
  
  muY.flat <- matrix(unlist(SIMP.fit$muY), ncol = length(SIMP.fit$muY[[1]]), byrow = T)
  muY.sd <- apply(muY.flat, 2, sd)
  muY.acf <- apply(muY.flat, 2, function(x) acf(x, lag.max = acf.lagmax))
  
  beta1C.flat <- matrix(unlist(SIMP.fit$beta1C), ncol = prod(dim(SIMP.fit$beta1C[[1]])), byrow = T)
  beta1C.sd <- apply(beta1C.flat, 2, sd)
  beta1C.acf <- apply(beta1C.flat, 2, function(x) acf(x, lag.max = acf.lagmax))                   
  
  
  beta1D.flat <- matrix(unlist(SIMP.fit$beta1D), ncol = prod(dim(SIMP.fit$beta1D[[1]])), byrow = T)
  beta1D.sd <- apply(beta1D.flat, 2, sd)
  beta1D.acf <- apply(beta1D.flat, 2, function(x) acf(x, lag.max = acf.lagmax)) 
  
  beta2.flat <- matrix(unlist(SIMP.fit$beta2), ncol = prod(dim(SIMP.fit$beta2[[1]])), byrow = T)
  beta2.sd <- apply(beta2.flat, 2, sd)
  beta2.acf <- apply(beta2.flat, 2, function(x) acf(x, lag.max = acf.lagmax)) 

  
  SigmaCD.flat <- matrix(unlist(SIMP.fit$SigmaCD), ncol = prod(dim(SIMP.fit$SigmaCD[[1]])), byrow = T)
  SigmaCD.sd <- apply(SigmaCD.flat, 2, sd)
  SigmaCD.acf <- apply(SigmaCD.flat, 2, function(x) acf(x, lag.max = acf.lagmax)) 
  
  SigmaYX.flat <- matrix(unlist(SIMP.fit$SigmaYX), ncol = prod(dim(SIMP.fit$SigmaYX[[1]])), byrow = T)
  SigmaYX.sd <- apply(SigmaYX.flat, 2, sd)
  SigmaYX.acf <- apply(SigmaYX.flat, 2, function(x) acf(x, lag.max = acf.lagmax)) 
  
  sd <- list(muX1C.sd = muX1C.sd,
             muY.sd = muY.sd,
             beta1C.sd = beta1C.sd,
             beta1D.sd = beta1D.sd,
             beta2.sd = beta2.sd,
             SigmaCD.sd = SigmaCD.sd,
             SigmaYX.sd = SigmaYX.sd)
  
  acf <-  list(muX1C.acf = muX1C.acf,
               muY.acf = muY.acf,
               beta1C.acf = beta1C.acf,
               beta1D.acf = beta1D.acf,
               beta2.acf = beta2.acf,
               SigmaCD.acf = SigmaCD.acf,
               SigmaYX.acf = SigmaYX.acf)

  list(sd=sd, acf=acf)
}
