#' Function to compute A from from
#' @param gamma The basis matrix gamma of the envelope space.
#'
#' @export
find_A_from_gamma <- function(gamma) {
  m <- ncol(gamma)
  G1 <- as.matrix(gamma[1:m, ])
  # check if G1 is invertible - else reorganize the predictors
  Y.order <- 1:nrow(gamma)
  switch.flag <- FALSE
  if (abs(det(G1)) < 1e-30) {
    gamma.t <- t(gamma)
    Y.order <- qr(gamma.t, tol = 1e-7)$pivot
    #X <- X[, X.order]
    gamma <- gamma[Y.order, ]
    switch.flag <- TRUE
  }
  if (sum(Y.order!=(1:nrow(gamma))) == 0){
    switch.flag <- FALSE
  }
  G1 <- as.matrix(gamma[1:m, ])
  G2 <- gamma[-(1:m), ]
  A <- G2 %*% solve(G1)
  list(A = A, switch.flag = switch.flag, Y.order = Y.order)
}

#' Function to compute gamma from A
#' @param A The matrix A
#'
#' @export
find_gammas_from_A <- function(A) {
  dims <- dim(A)
  u <- dims[2]
  r <- sum(dims)
  CA <- matrix(0, nrow = r, ncol = u)
  DA <- matrix(0, nrow = r, ncol = r-u)
  CA[(u+1):r, ] <- A
  CA[1:u, 1:u] <- diag(1, u)
  DA[1:u, ] <- -t(A)
  DA[-(1:u), ] <- diag(1, r-u)
  CAtCA <- crossprod(CA)
  DAtDA <- crossprod(DA)
  gamma <- CA %*% sqrtmatinv(CAtCA)
  gamma0 <- DA %*% sqrtmatinv(DAtDA)

  list(gamma = gamma, gamma0 = gamma0,
       CA = CA, CAtCA = CAtCA,
       DA = DA, DAtDA = DAtDA)
}
