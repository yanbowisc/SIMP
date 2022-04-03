#' Bura_Cook estimator for the estimation of rank(beta_{1C})
#'
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X1D Design matrix of the discrete part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X2 Design matrix of the nuisance predictors. Must have the same number of rows as \code{Y}.
#' @param Y Response matrix. Must have the same number of rows as \code{X}.
#' @param sig.level Significance level of the sequence of Chi-squared tests called in 
#' Bura-Cook estimator.
#' @examples
#' \dontrun{
#' library(SIMP)
#' library(Renvlp)  # Load Renvlp package only for wheatprotein dataset.
#' data(wheatprotein)
#' set.seed(1)
#' d <- Bura_Cook(X1C = wheatprotein[, 4:5], X1D = matrix(wheatprotein[, 8], ncol = 1), 
#'                 X2 = wheatprotein[, 6:7], Y = wheatprotein[, 1:3], sig.level = 0.05)
#' }
#' @export
Bura_Cook <- function(X1C, X1D, X2, Y, sig.level = 0.05){
  
  Y <- as.matrix(Y)
  samp.size <- nrow(Y)
  r <- ncol(Y)
  Y_bar <- colMeans(Y)
  Y_ctr <- sweep(Y, 2, Y_bar, "-")
  
  X1C <- as.matrix(X1C)
  pc <- ncol(X1C)
  X1C_bar <- colMeans(X1C)
  X1C_ctr <- sweep(X1C, 2, X1C_bar, "-")
  
  if (!is.null(X1D)){
    X1D <- as.matrix(X1D)
    pd <- ncol(X1D)
    X1D_bar <- colMeans(X1D)
    X1D_ctr <- sweep(X1D, 2, X1D_bar, "-")
  }else{
    pd <- 0
  }
  
  if (!is.null(X2)){
    X2 <- as.matrix(X2)
    p2 <- ncol(X2)
    X2_bar <- colMeans(X2)
    X2_ctr <- sweep(X2, 2, X2_bar, "-")
  }else{
    p2 <- 0
  }
  
  if ((pd > 0) & (p2 > 0)){
    ols.fit.tmp <- lm(Y_ctr ~ X1D_ctr + X2_ctr)
    ols.fit <- lm(ols.fit.tmp$residuals ~ X1C_ctr)
  }else if((pd > 0) & (p2 == 0)){
    ols.fit.tmp <- lm(Y_ctr ~ X1D_ctr)
    ols.fit <- lm(ols.fit.tmp$residuals ~ X1C_ctr)
  }else if((pd == 0) & (p2 > 0)){
    ols.fit.tmp <- lm(Y_ctr ~ X2_ctr)
    ols.fit <- lm(ols.fit.tmp$residuals ~ X1C_ctr)
  }else{
    ols.fit <- lm(Y_ctr ~ X1C_ctr)
  }
  bura_std <- sqrt((samp.size - pc - 1)/samp.size) * sqrtmat(cov(X1C_ctr) * (samp.size - 1)/samp.size) %*% ols.fit$coefficients[-1, ] %*% sqrtmatinv(cov(ols.fit$residuals) * (samp.size - 1)/samp.size)
  
  rank_max <- min(pc, r)
  Lambda <- numeric(rank_max)
  bura_p_value <- numeric(rank_max)
  bura_phi_square <- (svd(bura_std)$d)^2
  for (k in 0:(rank_max - 1)){
    Lambda[k+1] <- samp.size * sum(bura_phi_square[(k+1):rank_max])
    bura_p_value[k+1] <- pchisq(Lambda[k+1], (pc - k) * (r - k), lower.tail = F)
  }
  if (all( bura_p_value < sig.level)){
    rank <- rank_max
  }else{
    rank <- which(bura_p_value > sig.level)[1] - 1
  }
  
  rank
}
