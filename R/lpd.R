#' Function to compute the log full conditional density of A
#' @param A Parameter A from previous iteration.
#' @param X1C_ctr Centered X1C.
#' @param X1C_ctr_shifted A inner variable.
#' @param Yctr Centered Y.
#' @param etaC Updated etaC in this iteration.
#' @param R R from previous iteration.
#' @param dy dy.
#' @param pd pd.
#' @param X1D_ctr_etaD.t A inner variable.
#' @param X2_ctr_beta2 A inner variable.
#' @param Omega.inv Inverse of the updated Omega in this iteration.
#' @param Omega0.inv Inverse of the updated Omega0 in this iteration.
#' @param SigmaYX.inv Inverse of the updated SigmaYX in this iteration.
#' @param gamma_shifted A inner variable.
#' @param Lambda.half Inverse of Lambda.
#' @param A0 A0.
#' @param KA.half.inv The (-1/2)-th power of KA.
#' @param SigmaA.half.inv The (-1/2)-th power of SigmaA.
#'
#' @export
lpd_A <- function(A,
                  X1C_ctr,
                  X1C_ctr_shifted,
                  Yctr,
                  etaC,
                  R,
                  dy,
                  pd,
                  X1D_ctr_etaD.t,
                  X2_ctr_beta2,
                  Omega.inv,
                  Omega0.inv,
                  SigmaYX.inv,
                  gamma_shifted,
                  Lambda.half,
                  A0,
                  KA.half.inv,
                  SigmaA.half.inv) {

L_L0 <- find_gammas_from_A(A)
L <- L_L0$gamma
L0 <- L_L0$gamma0
SigmaCD.inv <- L %*% tcrossprod(Omega.inv, L) + L0 %*% tcrossprod(Omega0.inv, L0)
if (dy == 0){
  resi <- Yctr - X2_ctr_beta2
}else{
  resi <- Yctr - tcrossprod(tcrossprod(X1C_ctr%*%L, etaC) + X1D_ctr_etaD.t, R) - X2_ctr_beta2
}
term1 <- sum(SigmaCD.inv*crossprod(X1C_ctr_shifted))
term2 <- sum(SigmaYX.inv*crossprod(resi))
if (pd > 0){
  term3 <- sum(SigmaCD.inv*crossprod(Lambda.half%*%gamma_shifted))
}else{
  term3 <- 0
}
term4 <- sum((KA.half.inv %*% (A-A0) %*% SigmaA.half.inv)^2)

out <- - 0.5 * (term1 + term2 + term3 + term4)
attr(out, "L_L0") <- L_L0
out
}

#' Function to compute the log full conditional density of B
#' @param A Parameter B from previous iteration.
#' @param X1C_ctr Centered X1C.
#' @param Yctr Centered Y.
#' @param L Updated L.
#' @param dx dx
#' @param pd pd
#' @param etaC Updated etaC in this iteration.
#' @param X1D_ctr_etaD.t A inner variable.
#' @param X2_ctr_beta2 A inner variable.
#' @param Phi.inv Inverse of the updated Phi in this iteration.
#' @param Phi0.inv Inverse of the updated Phi0 in this iteration.
#' @param beta2.shifted A inner variable.
#' @param M.half Square root of M.
#' @param B0 B0
#' @param KB.half.inv The (-1/2)-th power of KB.
#' @param SigmaB.half.inv The (-1/2)-th power of SigmaB.
#'
#' @export
lpd_B <- function(A,
                  X1C_ctr,
                  Yctr,
                  L,
                  dx,
                  pd,
                  etaC,
                  X1D_ctr_etaD.t,
                  X2_ctr_beta2,
                  Phi.inv,
                  Phi0.inv,
                  beta2.shifted,
                  M.half,
                  B0,
                  KB.half.inv,
                  SigmaB.half.inv) {
  if (is.null(nrow(beta2.shifted))){
    p2 <- 0
  }else{
    p2 <- nrow(beta2.shifted)
  }

  R_R0 <- find_gammas_from_A(A)
  R <- R_R0$gamma
  R0 <- R_R0$gamma0
  SigmaYX.inv <- R %*% tcrossprod(Phi.inv, R) + R0 %*% tcrossprod(Phi0.inv, R0)
  if (dx > 0){
    resi <- Yctr - tcrossprod(tcrossprod(X1C_ctr%*%L, etaC) + X1D_ctr_etaD.t, R) - X2_ctr_beta2
  }else{
    if (pd > 0){
      resi <- Yctr - tcrossprod(X1D_ctr_etaD.t, R) - X2_ctr_beta2
    }else{
      resi <- Yctr - X2_ctr_beta2
    }

  }

  term1 <- sum(SigmaYX.inv*crossprod(resi))
  if (p2 > 0){
    term2 <- sum(SigmaYX.inv*crossprod(M.half%*%beta2.shifted))
  }else{
    term2 <- 0
  }

  term3 <- sum((KB.half.inv %*% (A-B0) %*% SigmaB.half.inv)^2)

  out <- - 0.5 * (term1 + term2 + term3)
  attr(out, "R_R0") <- R_R0
  out
}
