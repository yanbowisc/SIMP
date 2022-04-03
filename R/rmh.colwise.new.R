#' Function to implement the Metropolis step.
#' @param A_start A or B in the current iteration.
#' @param lpd_func log full conditional density function of A or B.
#' @param method Method for the metropolis. We have three choices: RW (default), HMC, NUTS. We recommend to
#' use RW method always if no other reasons.
#' @param tau The standard deviation for the proposal distribution used in the RW-Metropolis.
#' @param eps Initial value for estimating the stepsize in NUTS-metropolis.
#' @param eps_bar eps_bar parameter for NUTS-metropolis.
#' @param H The H_m parameter for NUTS-metropolis.
#' @param mu The mu parameter for NUTS-metropolis.
#' @param LL Number of leapfrog steps in the trajectory
#' @param ColSigma.A.used Covariance matrix for our matrix parameter A or B
#' @param alpha Tempering parameter
#' @param grad_lpd_func The function to calculate the gradients of 
#' the log full conditional posterior density of A or B.
#' @param M_adapt The M_adapt parameter used in the Algorithm 6 of NUTS.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param delta The delta parameter for NUTS-metropolis.
#' @param iter The iteration index.
#' @param samp.size Sample size.
#' @param ... Other parameters that may be useful.
#'
#' @export
rmh_colwise_new <- function(A_start, lpd_func, method = "RW", tau = NULL, eps, eps_bar, H, mu, LL = 10,
                            ColSigma.A.used, alpha, grad_lpd_func, M_adapt, M_diag, delta, iter, samp.size,
                            ...) {

  if (!method %in% c("RW", "HMC", "NUTS"))
    stop("method must be either RW, HMC or NUTS!!")

  dims <- dim(A_start)
  u <- dims[2] # stands for the envelope dim
  v <- sum(dims) # stands for the maximal envelope dim
  accpt.A <- numeric(u)

  A <- A_start
  lpd <- lpd_func(A = A, ...)

  if (iter == 1){
    eps <- eps_bar <- H <- mu <- matrix(1, dims[1], dims[2])
  }


  # if (autotune_size == "single") tau_curr <- tau
  if (method == "RW"){
    for(j in sample(1:u)) {
      tau_curr <- tau[, j]
      A_j_star <- A[, j] + rnorm(v-u, 0, tau_curr)
      A_star <- A
      A_star[, j] <- A_j_star

      lpd_star <- lpd_func(A = A_star, ...)

      if (log(runif(1)) < (lpd_star - lpd)) {
        # accept
        A <- A_star
        lpd <- lpd_star
        accpt.A[j] <- 1
      }}
  }else if(method == "HMC") {
    for(j in sample(1:u)) {
      tau_curr <- tau[, j]
      HMC_outp <- HMC_new(current_A = A,
                      j = j,
                      lpd_func = lpd_func,
                      epsilon = tau_curr,
                      LL = LL,
                      ColSigma.A.used = ColSigma.A.used,
                      alpha = alpha,
                      current_U = -lpd, 
                      samp.size = samp.size,
                      ...)
      A <- HMC_outp$A
      lpd <- HMC_outp$lpd
      accpt.A[j] <- HMC_outp$accept
    }
  }else if(method == "NUTS") {
    for(j in sample(1:u)) {

      NUTS_outp <- NUTS(A = A,
                       j = j,
                       iter = iter,
                       lpd_func = lpd_func,
                       grad_lpd_func = grad_lpd_func,
                       M_adapt = M_adapt,
                       M_diag = M_diag,
                       delta = delta,
                       max_treedepth = 6,
                       eps = eps[1, j],
                       eps_bar = eps_bar[1, j],
                       H = H[1, j],
                       mu = mu[1, j],
                       verbose = TRUE,
                       Delta_max = 1000,
                       ...)
      A <- NUTS_outp$A
      lpd <- NUTS_outp$lpd
      accpt.A[j] <- NUTS_outp$accept
      eps[, j] <- NUTS_outp$eps
      eps_bar[, j] <- NUTS_outp$eps_bar
      H[, j] <- NUTS_outp$H
      mu[, j] <- NUTS_outp$mu
    }

  }

  out <- list(A = A,
              lpd = lpd,
              accpt = accpt.A,
              tau = tau,
              eps = eps,
              eps_bar = eps_bar,
              H = H,
              mu = mu)

}

