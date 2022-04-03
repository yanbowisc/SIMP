#' Calculation of the posterior credible intervals for regression parameters
#'
#' @param SIMP.fit The MCMC output from \code{SIMP()} function.
#' @param combined Logical. Indicate whether posterior samples (after burn-in) from all chains 
#' should be combined to calculate the posterior credible interval. If there is only one chain,
#' then there is no difference by indicating this parameter to be TRUE or FALSE.
#' @param Chain.no If combined = FALSE, we need to indicate which chain should be used for the 
#' calculation of posterior credible interval. Chain.no must be an integer between 1 and the number of chains.
#' @param levels Significance level, default to be 0.95.
#' @examples
#' \dontrun{
#' library(SIMP)
#' library(Renvlp)
#' data(wheatprotein) # Load Renvlp package only for wheatprotein dataset.
#' set.seed(1)
#' X1C = wheatprotein[, 4:5]
#' X1D = as.matrix(wheatprotein[, 8], ncol = 1)
#' X2 = wheatprotein[, 6:7]
#' Y = wheatprotein[, 1:3]
#' MC_output <- SIMP(X1C = X1C, X1D = X1D, X2 = X2,
#'                   Y = Y, dx = 1, dy = 1, n.iter = 1e4)
#' CI(MC_output)
#' }
#' @export
CI <- function(SIMP.fit, combined = T, Chain.no = 1, levels = 0.95){
  
  pos.samp.size <- length(SIMP.fit$muX1C[[1]])
  n.chains <- length(SIMP.fit$muX1C)
  
  if (combined){
    muX1C.combined <- do.call(c, SIMP.fit$muX1C)
    names(muX1C.combined) <- NULL
  }else{
    muX1C.combined <- SIMP.fit$muX1C[[Chain.no]]
  }
  
  muX1C.flat <- matrix(unlist(muX1C.combined), ncol = length(muX1C.combined[[1]]), byrow = T)
  muX1C.lower <- apply(muX1C.flat, 2, function(x) quantile(x, probs = (1 - levels)/2))
  muX1C.upper <- apply(muX1C.flat, 2, function(x) quantile(x, probs = (1 + levels)/2))
  
  if (combined){
    muY.combined <- do.call(c, SIMP.fit$muY)
    names(muY.combined) <- NULL
  }else{
    muY.combined <- SIMP.fit$muY[[Chain.no]]
  }
  
  muY.flat <- matrix(unlist(muY.combined), ncol = length(muY.combined[[1]]), byrow = T)
  muY.lower <- apply(muY.flat, 2, function(x) quantile(x, probs = (1 - levels)/2))
  muY.upper <- apply(muY.flat, 2, function(x) quantile(x, probs = (1 + levels)/2))
  
  
  if (combined){
    beta1C.combined <- do.call(c, SIMP.fit$beta1C)
    names(beta1C.combined) <- NULL
  }else{
    beta1C.combined <- SIMP.fit$beta1C[[Chain.no]]
  }
  beta1C.flat <- matrix(unlist(beta1C.combined), ncol = prod(dim(beta1C.combined[[1]])), byrow = T)
  beta1C.lower <- matrix(apply(beta1C.flat, 2, function(x) quantile(x, probs = (1 - levels)/2)), ncol = ncol(beta1C.combined[[1]]))
  beta1C.upper <- matrix(apply(beta1C.flat, 2, function(x) quantile(x, probs = (1 + levels)/2)), ncol = ncol(beta1C.combined[[1]]))
  
  if (combined){
    beta1D.combined <- do.call(c, SIMP.fit$beta1D)
    names(beta1D.combined) <- NULL
  }else{
    beta1D.combined <- SIMP.fit$beta1D[[Chain.no]]
  }
  beta1D.flat <- matrix(unlist(beta1D.combined), ncol = prod(dim(beta1D.combined[[1]])), byrow = T)
  beta1D.lower <- matrix(apply(beta1D.flat, 2, function(x) quantile(x, probs = (1 - levels)/2)), ncol = ncol(beta1D.combined[[1]]))
  beta1D.upper <- matrix(apply(beta1D.flat, 2, function(x) quantile(x, probs = (1 + levels)/2)), ncol = ncol(beta1D.combined[[1]]))
  
  if (combined){
    beta2.combined <- do.call(c, SIMP.fit$beta2)
    names(beta2.combined) <- NULL
  }else{
    beta2.combined <- SIMP.fit$beta2[[Chain.no]]
  }
  beta2.flat <- matrix(unlist(beta2.combined), ncol = prod(dim(beta2.combined[[1]])), byrow = T)
  beta2.lower <- matrix(apply(beta2.flat, 2, function(x) quantile(x, probs = (1 - levels)/2)), ncol = ncol(beta2.combined[[1]]))
  beta2.upper <- matrix(apply(beta2.flat, 2, function(x) quantile(x, probs = (1 + levels)/2)), ncol = ncol(beta2.combined[[1]]))
  
  if (combined){
    SigmaCD.combined <- do.call(c, SIMP.fit$SigmaCD)
    names(SigmaCD.combined) <- NULL
  }else{
    SigmaCD.combined <- SIMP.fit$SigmaCD[[Chain.no]]
  }
  SigmaCD.flat <- matrix(unlist(SigmaCD.combined), ncol = prod(dim(SigmaCD.combined[[1]])), byrow = T)
  SigmaCD.lower <- matrix(apply(SigmaCD.flat, 2, function(x) quantile(x, probs = (1 - levels)/2)), ncol = ncol(SigmaCD.combined[[1]]))
  SigmaCD.upper <- matrix(apply(SigmaCD.flat, 2, function(x) quantile(x, probs = (1 + levels)/2)), ncol = ncol(SigmaCD.combined[[1]]))
  
  if (combined){
    SigmaYX.combined <- do.call(c, SIMP.fit$SigmaYX)
    names(SigmaYX.combined) <- NULL
  }else{
    SigmaYX.combined <- SIMP.fit$SigmaYX[[Chain.no]]
  }
  SigmaYX.flat <- matrix(unlist(SigmaYX.combined), ncol = prod(dim(SigmaYX.combined[[1]])), byrow = T)
  SigmaYX.lower <- matrix(apply(SigmaYX.flat, 2, function(x) quantile(x, probs = (1 - levels)/2)), ncol = ncol(SigmaYX.combined[[1]]))
  SigmaYX.upper <- matrix(apply(SigmaYX.flat, 2, function(x) quantile(x, probs = (1 + levels)/2)), ncol = ncol(SigmaYX.combined[[1]]))
  
  CI <- list(muX1C.lower = muX1C.lower,
             muX1C.upper = muX1C.upper,
             muY.lower = muY.lower,
             muY.upper = muY.upper,
             beta1C.lower = beta1C.lower,
             beta1C.upper = beta1C.upper,
             beta1D.lower = beta1D.lower,
             beta1D.upper = beta1D.upper,
             beta2.lower = beta2.lower,
             beta2.upper = beta2.upper,
             SigmaCD.lower = SigmaCD.lower,
             SigmaCD.upper = SigmaCD.upper,
             SigmaYX.lower = SigmaYX.lower,
             SigmaYX.upper = SigmaYX.upper)
  
  CI
}
