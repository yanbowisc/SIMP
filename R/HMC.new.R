#' Hamiltonian monte-carlo
#' @import madness
#' @param current_A Current iterate of our matrix parameter A or B
#' @param j Index for which column to update
#' @param lpd_func Log full conditional posterior density
#' @param epsilon Step-size
#' @param LL Number of leapfrog steps in the trajectory
#' @param ColSigma.A.used Covariance matrix for our matrix parameter A or B
#' @param alpha Tempering parameter
#' @param current_U Energy function of the current iterate of our matrix parameter A or B
#' @param samp.size Sample size
#' @param ... Other parameters needed for HMC updating.
#'
#' @export
HMC_new <- function(current_A, j, lpd_func, epsilon, LL, ColSigma.A.used = ColSigma.A.used, alpha, current_U, samp.size,
                    ...){
  dims <- dim(current_A)
  A <- current_A
  #AA_j <- rnorm(dims[1], 0, 1)
  AA_j <- rmvnorm(n = 1, mu = rep(0, dims[1]), sigma = solve_chol(ColSigma.A.used))
  current_AA_j <- AA_j

  grad_A <- -numderiv(f = lpd_func,
                      x = A,
                      eps = 1e-8,
                      type = "central", ...)/samp.size

  grad_A_j <- grad_A[((j-1)*dims[1] + 1) : (j*dims[1])]

  #epsilon <- runif(1, epsilon_lower, epsilon_upper)
  epsilon <- runif(dims[1], 0.9*epsilon, 1.1*epsilon)
  #  Make a half step for momentum at the beginning
  AA_j <- AA_j * sqrt(alpha)
  AA_j <- AA_j - epsilon * grad_A_j/2

  # Alternate full steps for positions and momentum
  for (l in 1:(LL/2)){
    #epsilon <- runif(1, epsilon_lower, epsilon_upper)

    # Make a full step for the position
    A[, j] <- A[, j] + epsilon * ColSigma.A.used %*% AA_j

    grad_A <- -numderiv(f = lpd_func,
                        x = A,
                        eps = 1e-8,
                        type = "central", ...)/samp.size
    grad_A_j <- grad_A[((j-1)*dims[1] + 1) : (j*dims[1])]

    # Make a full step for the momentum, except at end of trajectory
    if (l != (LL/2)) {
      AA_j <- AA_j - epsilon * grad_A_j
      AA_j <- AA_j*alpha}
  }

  # Make a half step for momentum at the end.
  AA_j <- AA_j - epsilon * grad_A_j

  for (l in (LL/2 + 1):LL){
    #epsilon <- runif(1, epsilon_lower, epsilon_upper)

    # Make a full step for the position
    A[, j] <- A[, j] + epsilon * ColSigma.A.used %*% AA_j

    grad_A <- -numderiv(f = lpd_func,
                        x = A,
                        eps = 1e-8,
                        type = "central", ...)/samp.size
    grad_A_j <- grad_A[((j-1)*dims[1] + 1) : (j*dims[1])]

    # Make a full step for the momentum, except at end of trajectory
    if (l != LL/2) {
      AA_j <- AA_j - epsilon * grad_A_j
      AA_j <- AA_j/alpha
      }
  }

  AA_j <- AA_j - epsilon * grad_A_j
  AA_j <- AA_j/sqrt(alpha)

  # Negate momentum at end of trajectory to make the proposal symmetric
  AA_j <- - AA_j

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_K <- crossprod(current_AA_j, ColSigma.A.used) %*% current_AA_j/2
  proposed_K <- crossprod(AA_j, ColSigma.A.used) %*% AA_j/2


  proposed_U <- -lpd_func(A, ...)

  # Accept or reject the state at the end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (log(runif(1)) < (current_U - proposed_U + current_K - proposed_K ))
  {
    # accept
    accept <- 1
    lpd <- - proposed_U
  }else{
    # reject
    A <- current_A
    accept <- 0
    lpd <- - current_U
  }
  return(list(A = A,
              accept = accept,
              lpd = lpd))

}
