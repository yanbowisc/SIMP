#' The function to generate a Inverse-Wishart sample
#' @param dim Dimension.
#' @param Phi Scale matrix parameter.
#' @param nu Degrees of freedom.
#'
#' @export
rinvwish <- function(dim, Phi, nu) {
  p <- dim
  Phi.inv <- (solve_chol(Phi))
  Sigma.inv <- rWishart(1, nu, Phi.inv)[, , 1]
  (solve_chol(Sigma.inv))
}

#' The function to generate vector Normal samples
#' @param n Sample size.
#'
#' @param mu Mean.
#' @param sigma Covariance matrix.
#' @param dim Dimension.
#'
#' @export
rmvnorm <- function(n = 1, mu, sigma, dim = NULL) {
  if (is.null(dim)) {
    dim <- ncol(sigma)
  }
  p <- dim
  # eig <- eigen(Sigma)
  # Sigma.half <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Sigma.half <- t(chol(sigma))
  z <- rnorm(p)
  as.vector(mu + as.vector(Sigma.half %*% z))
}

#' An efficient function to calculate the inverse of a matrix
#' @param mat The matrix that you want to calculate the inverse.
#' @export
solve_chol <- function(mat) chol2inv(chol(mat))

#' The function to calculate the square root of a matrix
#' @param mat The matrix that you want to calculate the square root.
#' @export
sqrtmat <- function(mat) {
  e <- eigen(mat, symmetric = TRUE)
  if(length(e$values) == 1)
    sqrt(mat)
  else
    e$vectors %*% diag(c(sqrt(e$values))) %*% t(e$vectors)
  # t(chol(mat))
}

#' The function to calculate the -(1/2)-th power of a matrix
#' @param mat The matrix that you want to calculate the -(1/2)-th power.
#' @export
sqrtmatinv <- function(mat) {
  e <- eigen(mat, symmetric = TRUE)
  if(length(e$values) == 1)
    1/sqrt(mat)
  else
    e$vectors %*% diag(1/c(sqrt(e$values))) %*% t(e$vectors)
  # t(chol(solve_chol(mat)))
}

#' The function to generate a Matrix Normal sample.
#' @param M Location parameter.
#'
#' @param V1 Scale parameter for the row.
#' @param V2 Scale parameter for the column.
#'
#' @export
rMatrixNormal<-function(M, V1, V2){
  #M is the mean matrix, and the V1, V2 are symmetric positive definite matrices

  alldims <- dim(M)
  n.rows <- alldims[1]
  n.cols <- alldims[2]
  z <- matrix(rnorm(n.rows * n.cols), nrow = n.rows, ncol = n.cols)

  A <- t(chol(V1))  # V1 = A %*% t(A)
  B <- chol(V2)  # V2 = t(B) %*% B
  M + A %*% z %*% B
}

#' Function for the progressbar
#' @importFrom matrixStats rowVars
#' @param ... Settings for the progress bar.
#' @export
tkProgressBar <- function(...) utils::txtProgressBar(..., style = 3)
setTkProgressBar <- utils::setTxtProgressBar
Matrix <- Matrix::Matrix
