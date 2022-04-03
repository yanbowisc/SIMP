#' Hyper-parameters specification function
#' @description either generating a set of hyper-parameters for un-informative priors, or 
#' passing a set of user-specified hyper-parameters to \code{SIMP()} function.
#' @param pc Dimension of X1C, the continuous part of the predictors of interest.
#' @param pd Dimension of X1D, the discrete part of the predictors of interest.
#' @param p2 Dimension of X2, predictors of not main interest.
#' @param r Dimension of response Y.
#' @param dx Partial predictor envelope dimension of SIMP.
#' @param dy Partial response envelope dimension of SIMP.
#' @param ... Input for the user-specified hyper-parameters.
#' @examples
#' \dontrun{
#' library(SIMP)
#' pd = 8; pd = 2; p2 = 2; r = 8; dx = 1; dy = 1
#' hp <- get_hyperparams_init(pc, pd, p2, r, dx, dy)
#' }
#' @export
get_hyperparams_init <- function(pc, pd, p2, r, dx, dy, ...){

  dots <- list(...)
  hp <- dots$hyper_params

  if (is.null(hp)) {
    hp = list()
  }

  if (is.null(hp$e)) { hp$e <- matrix(0, p2, r) }

  if (is.null(hp$g)) { hp$g <- matrix(0, pd, pc) }

  if (is.null(hp$h)) { hp$h <- matrix(0, dy, dx) }

  if (is.null(hp$f)) { hp$f <- matrix(0, dy, pd) }

  if (is.null(hp$M)) { hp$M <- 1e-6 * diag(1, nrow = p2) }

  if (is.null(hp$Lambda)) { hp$Lambda <- 1e-6 * diag(1, nrow = pd) }

  if (is.null(hp$E)) { hp$E <- 1e-6 * diag(1, nrow = dx) }

  if (is.null(hp$Q)) { hp$Q <- 1e-6 * diag(1, nrow = pd) }

  if (is.null(hp$PsiX)) { hp$PsiX <- 1e-6 * dx * diag(1, nrow = dx) }

  if (is.null(hp$wX)) { hp$wX <- dx }

  if (is.null(hp$PsiX0)) { hp$PsiX0 <- 1e-6 * (pc - dx) * diag(1, nrow = pc - dx) }

  if (is.null(hp$wX0)) { hp$wX0 <- pc - dx }

  if (is.null(hp$PsiY)) { hp$PsiY <- 1e-6 * dy * diag(1, nrow = dy) }

  if (is.null(hp$wY)) { hp$wY <- dy }

  if (is.null(hp$PsiY0)) { hp$PsiY0 <- 1e-6 * ( r - dy) * diag(1, nrow = r - dy) }

  if (is.null(hp$wY0)) { hp$wY0 <- r - dy }

  if (dx > 0 & dx < pc) {

    if (is.null(hp$A0)) { hp$A0 <- matrix(0, pc - dx, dx) }

    if (is.null(hp$KA)) { hp$KA <- 1e6 * diag(1, nrow = pc - dx) }

    hp$KA.inv <- solve_chol(hp$KA)
    hp$KA.half.inv <- sqrtmatinv(hp$KA)

    if (is.null(hp$SigmaA)) { hp$SigmaA <- 1e6 * diag(1, nrow = dx) }

    hp$SigmaA.inv <- solve_chol(hp$SigmaA)
    hp$SigmaA.half.inv <- sqrtmatinv(hp$SigmaA)

  }


  if (dy > 0 & dy < r) {

    if (is.null(hp$B0)) { hp$B0 <- matrix(0, r - dy, dy) }

    if (is.null(hp$KB)) { hp$KB <- 1e6 * diag(1, nrow = r - dy) }

    hp$KB.inv <- solve_chol(hp$KB)
    hp$KB.half.inv <- sqrtmatinv(hp$KB)

    if (is.null(hp$SigmaB)) { hp$SigmaB <- 1e6 * diag(1, nrow = dy) }

    hp$SigmaB.inv <- solve_chol(hp$SigmaB)
    hp$SigmaB.half.inv <- sqrtmatinv(hp$SigmaB)

  }

  hp

}
