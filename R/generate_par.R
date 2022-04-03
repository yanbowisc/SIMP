#' Parameter generating function in Simulation section
#'
#' @param r Dimension of response Y.
#' @param pc Dimension of X1C, the continuous part of the predictors of interest.
#' @param pd Dimension of X1D, the discrete part of the predictors of interest.
#' @param p2 Dimension of X2, predictors of not main interest.
#' @param dx Partial predictor envelope dimension of SIMP.
#' @param dy Partial response envelope dimension of SIMP.
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
#' }
#' @export
generate_par <- function(r, pc, pd, p2, dx, dy) {

  A_exists <- (dx > 0) & (dx < pc)
  B_exists <- (dy > 0) & (dy < r)
  muY.tru <- runif(r, 0, 10)
  muX1C.tru <- runif(pc, 0, 10)
  if(pd > 0){
    gamma.tru <- matrix(runif(pd * pc, min = -2, max = 2), nrow = pd, ncol = pc)
  }else{
    gamma.tru <- 0
  }

  if (p2 > 0){
    beta2.tru <- matrix(runif(p2 * r, min = -2, max = 2), nrow = p2, ncol = r)
  }else{
    beta2.tru <- 0
  }

  if (dx * dy > 0){
    etaC.tru <- matrix(runif(dy * dx, min = -2, max = 2), nrow = dy, ncol = dx)
  }else{
    etaC.tru <- 0
  }

  if (dy * pd > 0){
    etaD.tru <- matrix(runif(dy * pd, min = -2, max = 2), nrow = dy, ncol = pd)
  }else{
    etaD.tru <- 0
  }

  if (A_exists){
    A <- A.tru <- matrix(runif((pc - dx) * dx, min = -1, max = 1), nrow = pc - dx, ncol = dx)
    L_L0.tru <- find_gammas_from_A(A)
    L.tru <- L_L0.tru$gamma
    L0.tru <- L_L0.tru$gamma0
    # Omega.tru <- diag(sort(runif(dx, 10, 20), decreasing = TRUE), ncol = dx, nrow = dx)
    # Omega0.tru <- diag(sort(runif(pc - dx, 5, 10), decreasing = TRUE), ncol = pc - dx, nrow = pc - dx)
    Omega.tru <- diag(10, nrow = dx, ncol = dx)
    Omega0.tru <- diag(1, nrow = pc - dx, ncol = pc - dx)
    SigmaCD1 <- L.tru %*% Omega.tru %*% t(L.tru)
    SigmaCD2 <- L0.tru %*% Omega0.tru %*% t(L0.tru)
  }else if (dx == 0){
    A <- A.tru <- 0
    L.tru <- 0
    L0.tru <- diag(1, pc)
    Omega.tru <- 0
    Omega0.tru <- diag(1, nrow = pc, ncol = pc)
    SigmaCD1 <- matrix(0, pc, pc)
    SigmaCD2 <- Omega0.tru
  }else if (dx == pc){
    A <- A.tru <- 0
    L.tru <- diag(1, pc)
    L0.tru <- 0
    Omega.tru <- diag(10, nrow = pc, ncol = pc)
    Omega0.tru <- 0
    SigmaCD1 <- Omega.tru
    SigmaCD2 <- matrix(0, pc, pc)
  }

  if (B_exists){
    B <- B.tru <- matrix(runif((r - dy) * dy, min = -1, max = 1), nrow = r - dy, ncol = dy)
    R_R0.tru <- find_gammas_from_A(B)
    R.tru <- R_R0.tru$gamma
    R0.tru <- R_R0.tru$gamma0
    # Phi.tru <- diag(sort(runif(dy, 0, 1), decreasing = TRUE), ncol = dy, nrow = dy)
    # Phi0.tru <- diag(sort(runif(r - dy, 5, 10), decreasing = TRUE), ncol = r - dy, nrow = r - dy)
    Phi.tru <- diag(1, nrow = dy, ncol = dy)
    Phi0.tru <- diag(5, nrow = r - dy, ncol = r - dy)
    SigmaYX1 <- R.tru %*% Phi.tru %*% t(R.tru)
    SigmaYX2 <- R0.tru %*% Phi0.tru %*% t(R0.tru)
    if (pd > 0){
      beta1D.tru <- t(R.tru %*% etaD.tru)
    }else{
      beta1D.tru <- 0
    }

  }else if (dy == 0){
    B <- B.tru <- 0
    R.tru <- 0
    R0.tru <- diag(1, r)
    Phi.tru <- 0
    Phi0.tru <- diag(5, nrow = r, ncol = r)
    SigmaYX1 <- matrix(0, r, r)
    SigmaYX2 <- Phi0.tru
    if (pd > 0){
      beta1D.tru <- matrix(0, nrow = pd, ncol = r)
    }else{
      beta1D.tru <- 0
    }

  }else if (dy == r){
    B <- B.tru <- 0
    R.tru <- diag(1, r)
    R0.tru <- 0
    Phi.tru <- diag(1, nrow = r, ncol = r)
    Phi0.tru <- 0
    SigmaYX1 <- Phi.tru
    SigmaYX2 <- matrix(0, r, r)
    if (pd > 0){
      beta1D.tru <- t(R.tru %*% etaD.tru)
    }else{
      beta1D.tru <- 0
    }
  }

  if ((dx > 0)&(dy > 0)){
    beta1C.tru <- L.tru %*% t(R.tru %*% etaC.tru)
  }else{
    beta1C.tru <- matrix(0, nrow = pc, ncol = r)
  }


  SigmaCD.tru <- SigmaCD1 + SigmaCD2
  SigmaYX.tru <- SigmaYX1 + SigmaYX2



  list(muX1C.tru = muX1C.tru,
       muY.tru = muY.tru,
       beta1C.tru = beta1C.tru,
       beta1D.tru = beta1D.tru,
       beta2.tru = beta2.tru,
       gamma.tru = gamma.tru,
       Omega.tru = Omega.tru,
       Omega0.tru = Omega0.tru,
       Phi.tru = Phi.tru,
       Phi0.tru = Phi0.tru,
       A.tru = A.tru,
       L.tru = L.tru,
       L0.tru = L0.tru,
       B.tru = B.tru,
       R.tru = R.tru,
       R0.tru = R0.tru,
       etaC.tru = etaC.tru,
       etaD.tru = etaD.tru,
       SigmaCD.tru = SigmaCD.tru,
       SigmaYX.tru = SigmaYX.tru)
}
