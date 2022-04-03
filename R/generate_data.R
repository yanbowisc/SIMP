#' Data generating function according to SIMP
#'
#' @param muX1C.tru a vector of length p_C. The true value of mu1C.
#' @param muY.tru a vector of length r. The true value of muY.
#' @param beta1C.tru a p_C by r matrix. The true value of beta1C.
#' @param beta1D.tru a p_D by r matrix. The true value of beta1D.
#' @param beta2.tru a p_2 by r matrix. The true value of beta2.
#' @param gamma.tru a p_D by p_C matrix. The true value of gamma.
#' @param SigmaCD.tru a p_C by p_C matrix. The true value of SigmaCD.
#' @param SigmaYX.tru a r by r matrix. The true value of SigmaYX.
#' @param K a positive integer. The K in the discrete uniform{0,1,...,K-1} for the generation of X1D.
#' @param mu2 a vector of length p_2. The true value of mean of Normal distribution for the generation of X2.
#' @param SigmaX2 a p_2 by p_2 matrix. The true value of Covariance matrix of Normal distribution for the generation of X2.
#' @param n Sample size.
#' @param r Dimension of response Y.
#' @param pc Dimension of X1C, the continuous part of the predictors of interest.
#' @param pd Dimension of X1D, the discrete part of the predictors of interest.
#' @param p2 Dimension of X2, predictors of not main interest.
#' @param ... Other parameters needed
#' @examples
#' \dontrun{
#' r <- 8
#' pc <- 8
#' pd <- 2
#' p2 <- 2
#' p = pc + pd + p2
#' K <- 3
#' mu2 <- c(2, 5)
#' dx.tru <- 6
#' dy.tru <- 2
#' n <- 300
#' set.seed(2)

#' if (p2 > 0){
#'   SigmaX2 <- rinvwish(p2, diag(1, p2), p2)
#' }else{
#'   SigmaX2 <- 0
#' }
#' all_pars <- generate_par(r, pc, pd, p2, dx.tru, dy.tru)
#'  dat <- do.call(generate_data, c(all_pars, 
#'  list(K = K, mu2 = mu2, SigmaX2 = SigmaX2, n = n, r = r, pc = pc, pd = pd, p2 = p2)))
#' }
#' @export
generate_data <- function(muX1C.tru, muY.tru,
                          beta1C.tru, beta1D.tru, beta2.tru,
                          gamma.tru,
                          SigmaCD.tru, SigmaYX.tru,
                          K, mu2, SigmaX2,
                          n, r, pc, pd, p2, ...)
{
  if (pd * p2 >0){
    # Generate X1D from Uniform dist.
    X1D <- matrix( floor(runif(n*pd, min = 0, max = K)), nrow = n, ncol = pd)
    X1D_ctr <- scale(X1D, center = TRUE, scale = FALSE)

    # Generate X1C from conditional normal on X1D
    epsCD <- matrix(rnorm(n*pc), n, pc) %*% sqrtmat(SigmaCD.tru)
    X1C <- tcrossprod(rep(1, n), muX1C.tru) + X1D_ctr %*% gamma.tru + epsCD
    X1C_ctr <- scale(X1C, center = TRUE, scale = FALSE)

    # Generate X2 from Normal
    X2 <- tcrossprod(rep(1, n), mu2) + matrix(rnorm(n*p2), n, p2) %*% sqrtmat(SigmaX2)
    X2_ctr <- scale(X2, center = TRUE, scale = FALSE)

    # Generate Y
    epsYX <- matrix(rnorm(n*r), n, r) %*% sqrtmat(SigmaYX.tru)
    Y <- tcrossprod(rep(1, n), muY.tru) + X1C_ctr%*%beta1C.tru + X1D_ctr%*%beta1D.tru + X2_ctr%*%beta2.tru + epsYX

    X <- cbind( X1C, X1D, X2)
  }else if ((pd == 0)&(p2 > 0)){
    # No X1D existed
    X1D <- NULL
    # Generate X1C from conditional normal on X1D
    epsCD <- matrix(rnorm(n*pc), n, pc) %*% sqrtmat(SigmaCD.tru)
    X1C <- tcrossprod(rep(1, n), muX1C.tru) + epsCD
    X1C_ctr <- scale(X1C, center = TRUE, scale = FALSE)

    # Generate X2 from Normal
    X2 <- tcrossprod(rep(1, n), mu2) + matrix(rnorm(n*p2), n, p2) %*% sqrtmat(SigmaX2)
    X2_ctr <- scale(X2, center = TRUE, scale = FALSE)

    # Generate Y
    epsYX <- matrix(rnorm(n*r), n, r) %*% sqrtmat(SigmaYX.tru)
    Y <- tcrossprod(rep(1, n), muY.tru) + X1C_ctr%*%beta1C.tru + X2_ctr%*%beta2.tru + epsYX

    X <- cbind( X1C, X2)
  }else if ((p2 == 0)&(pd > 0)){
    # Generate X1D from Uniform dist.
    X1D <- matrix( floor(runif(n*pd, min = 0, max = K)), nrow = n, ncol = pd)
    X1D_ctr <- scale(X1D, center = TRUE, scale = FALSE)

    # Generate X1C from conditional normal on X1D
    epsCD <- matrix(rnorm(n*pc), n, pc) %*% sqrtmat(SigmaCD.tru)
    X1C <- tcrossprod(rep(1, n), muX1C.tru) + X1D_ctr %*% gamma.tru + epsCD
    X1C_ctr <- scale(X1C, center = TRUE, scale = FALSE)

    # No X2 existed
    X2 <- NULL

    # Generate Y
    epsYX <- matrix(rnorm(n*r), n, r) %*% sqrtmat(SigmaYX.tru)
    Y <- tcrossprod(rep(1, n), muY.tru) + X1C_ctr%*%beta1C.tru + X1D_ctr%*%beta1D.tru + epsYX

    X <- cbind( X1C, X1D)
  }else if ((p2 == 0)&(pd == 0)){
    # No X1D existed.
    X1D <- NULL

    # Generate X1C from conditional normal on X1D
    epsCD <- matrix(rnorm(n*pc), n, pc) %*% sqrtmat(SigmaCD.tru)
    X1C <- tcrossprod(rep(1, n), muX1C.tru) + epsCD
    X1C_ctr <- scale(X1C, center = TRUE, scale = FALSE)

    # No X2 existed
    X2 <- NULL

    # Generate Y
    epsYX <- matrix(rnorm(n*r), n, r) %*% sqrtmat(SigmaYX.tru)
    Y <- tcrossprod(rep(1, n), muY.tru) + X1C_ctr%*%beta1C.tru + epsYX

    X <- X1C
  }

  list(X1C = X1C, X1D = X1D, X2 = X2, X = X, Y = Y)
}
