#' The function to calculate the gradients of the log full conditional posterior density of A or B.
#' @param A Parameter A or B
#'
#' @param j The index for which column of A or B to be updated.
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param ... Other parameters that may be useful.
#'
#' @export
grad_lpd_func <- function(A, j, lpd_func, ...){
  dims <- dim(A)
  grad_A <- numderiv(f = lpd_func,
                     x = A,
                     eps = 1e-8,
                     type = "central", ...)
  grad_A_j <- grad_A[((j-1)*dims[1] + 1) : (j*dims[1])]
  grad_A_j
}

#' The function to implement the leapfrop step
#' @param A Parameter A or B.
#'
#' @param aa The momentum vector of the chosen column of A or B.
#' @param eps Stepsize.
#' @param grad_lpd_func The function to calculate the gradients of 
#' the log full conditional posterior density of A or B.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param j The index for which column of A or B to be updated.
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param ... Other parameters that may be useful.
#'
#' @export
leapfrog_step = function(A, aa, eps, grad_lpd_func, M_diag, j, lpd_func, ...){
  A_tilde <- A
  aa_tilde <- aa + 0.5 * eps * grad_lpd_func(A, j, lpd_func, ...)
  A_tilde[, j] <- A[, j] + eps * M_diag %*% aa_tilde
  aa_tilde <- aa_tilde + 0.5 * eps * grad_lpd_func(A_tilde, j, lpd_func, ...)
  list(A = A_tilde, aa = aa_tilde)
}

#' The log joint full conditional density of the j-th column of A (or B) and its momentum vector.
#'
#' @param A Parameter A or B.
#' @param aa Momentum vector of the j-th column of A (or B).
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param ... Other parameters that may be useful.
#' @export
#'
joint_log_density = function(A, aa, lpd_func, M_diag, ...){
  as.numeric(lpd_func(A, ...)) - 0.5*crossprod(aa, M_diag%*%aa)
}

#' The function to find the stepsize
#' @param A Parameter A or B.
#' @param j The index for which column of A or B to be updated.
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param grad_lpd_func The function to calculate the gradients of 
#' the log full conditional posterior density of A or B.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param eps The initial value for the stepsize
#' @param verbose Whether to show the progress messages.
#' @param ... Other parameters that may be useful.
#'
#' @export
find_reasonable_epsilon = function(A, j, lpd_func, grad_lpd_func, M_diag, eps = 1, verbose = TRUE, ...){
  dims <- dim(A)
  aa <- rmvnorm(dim = dims[1], mu = 0, sigma = solve_chol(M_diag))
  proposed <- leapfrog_step(A, aa, eps, grad_lpd_func, M_diag, j, lpd_func, ...)
  log_ratio <- joint_log_density(proposed$A, proposed$aa, lpd_func, M_diag, ...) - joint_log_density(A, aa, lpd_func, M_diag, ...)
  alpha <- ifelse(exp(log_ratio) > 0.5, 1, -1)
  if(is.nan(alpha)) alpha <- -1
  count <- 1
  while(is.nan(log_ratio) || alpha * log_ratio > (-alpha)*log(2)){
    eps <- c(2^alpha * eps)
    proposed <- leapfrog_step(A, aa, eps, grad_lpd_func, M_diag, j, lpd_func, ...)
    log_ratio <- joint_log_density(proposed$A, proposed$aa, lpd_func, M_diag, ...) - joint_log_density(A, aa, lpd_func, M_diag, ...)
    count <- count + 1
    if(count > 100) {
      # stop("Could not find reasonable epsilon in 100 iterations!")
      break
    }
  }
  if(verbose) message("Reasonable epsilon = ", eps, " found after ", count-1, " steps")
  eps
}

#' Conditions checking function used in NUTS.
#' @param s An output from build_tree()
#'
#' @param A_plus An output from build_tree()
#' @param A_minus An output from build_tree()
#' @param aa_plus An output from build_tree()
#' @param aa_minus An output from build_tree()
#' @param j The index for which column of A or B to be updated.
#'
#' @export
check_NUTS = function(s, A_plus, A_minus, aa_plus, aa_minus, j){
  if(is.na(s)) return(0)
  condition1 <- crossprod(A_plus[, j] - A_minus[, j], aa_minus) >= 0
  condition2 <- crossprod(A_plus[, j] - A_minus[, j], aa_plus) >= 0
  s && condition1 && condition2
}

#' Tree building function for NUTS
#' @param A Parameter A or B.
#' @param j The index for which column of A or B to be updated.
#' @param aa The momentum vector of the j-th column of A or B.
#' @param u a random number used in NUTS.
#' @param v a random number from discrete Unif{-1, +1} used in NUTS.
#' @param jj The index for the layer of the tree.
#' @param eps Stepsize.
#' @param A00 A or B from the previous iteration.
#' @param aa0 The momentum vector from the previous iteration.
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param grad_lpd_func The function to calculate the gradients of 
#' the log full conditional posterior density of A or B.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param Delta_max The Delta_max parameter used in the Algorithm 6 of NUTS.
#' @param ... Other parameters that may be useful.
#'
#' @export
build_tree = function(A, j, aa, u, v, jj, eps, A00, aa0, lpd_func, grad_lpd_func, M_diag, Delta_max = 1000, ...){
  if(jj == 0){
    proposed <- leapfrog_step(A, aa, v*eps, grad_lpd_func, M_diag, j, lpd_func, ...)
    A <- proposed$A
    aa <- proposed$aa
    log_prob <- joint_log_density(A, aa, lpd_func, M_diag, ...)
    log_prob0 <- joint_log_density(A00, aa0, lpd_func, M_diag, ...)
    n <- (log(u) <= log_prob)
    s <- (log(u) < Delta_max + log_prob)
    alpha <- min(1, exp(log_prob - log_prob0))
    if(is.nan(alpha)) stop()
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(A_minus=A, A_plus=A, A=A, aa_minus=aa,
                aa_plus=aa, s=s, n=n, alpha=alpha, n_alpha=1))
  } else{
    obj0 <- build_tree(A=A, j=j, aa=aa, u=u, v=v, jj=jj-1, eps=eps, A00=A00, aa0=aa0, lpd_func=lpd_func, grad_lpd_func=grad_lpd_func, M_diag=M_diag, Delta_max=Delta_max, ...)
    A_minus <- obj0$A_minus
    aa_minus <- obj0$aa_minus
    A_plus <- obj0$A_plus
    aa_plus <- obj0$aa_plus
    A <- obj0$A
    if(obj0$s == 1){
      if(v == -1){
        obj1 <- build_tree(A=obj0$A_minus, j=j, aa=obj0$aa_minus, u=u, v=v, jj=jj-1, eps=eps, A00=A00, aa0=aa0, lpd_func=lpd_func, grad_lpd_func=grad_lpd_func, M_diag=M_diag, Delta_max=Delta_max, ...)
        A_minus <- obj1$A_minus
        aa_minus <- obj1$aa_minus
      } else{
        obj1 <- build_tree(A=obj0$A_plus, j=j, aa=obj0$aa_plus, u=u, v=v, jj=jj-1, eps=eps, A00=A00, aa0=aa0, lpd_func=lpd_func, grad_lpd_func=grad_lpd_func, M_diag=M_diag, Delta_max=Delta_max, ...)
        A_plus <- obj1$A_plus
        aa_plus <- obj1$aa_plus
      }
      n <- obj0$n + obj1$n
      if(n != 0){
        prob <- obj1$n / n
        if(runif(1) < prob){
          A <- obj1$A
        }
      }
      s <- check_NUTS(obj1$s, A_plus, A_minus, aa_plus, aa_minus, j)
      alpha <- obj0$alpha + obj1$alpha
      n_alpha <- obj0$n_alpha + obj1$n_alpha

    } else{
      n <- obj0$n
      s <- obj0$s
      alpha <- obj0$alpha
      n_alpha <- obj0$n_alpha
    }
    if(is.na(s) || is.nan(s)){
      s <- 0
    }
    if(is.na(n) || is.nan(n)){
      n <- 0
    }
    return(list(A_minus=A_minus, A_plus=A_plus, A=A,
                aa_minus=aa_minus, aa_plus=aa_plus, s=s, n=n, alpha=alpha, n_alpha=n_alpha))
  }
}

#' The function to implement The No-U-Turn Sampler (NUTS)
#' @param A Parameter A or B.
#' @param j The index for which column of A or B to be updated.
#' @param iter The iteration index.
#' @param lpd_func The log full conditional posterior density of A or B.
#' @param grad_lpd_func The function to calculate the gradients of 
#' the log full conditional posterior density of A or B.
#' @param M_adapt The M_adapt parameter used in the Algorithm 6 of NUTS.
#' @param M_diag Covariance matrix for our matrix parameter A or B.
#' @param delta The delta parameter used in the Algorithm 6 of NUTS.
#' @param max_treedepth The maximal tree depth.
#' @param eps The initial value for the stepsize tuning for find_reasonable_epsilon().
#' @param eps_bar The eps_bar parameter used in the Algorithm 6 of NUTS.
#' @param H The H_m in the Algorithm 6 of NUTS.
#' @param mu The mu in the Algorithm 6 of NUTS.
#' @param verbose Whether to show the progress messages for find_reasonable_epsilon().
#' @param Delta_max The Delta_max parameter used in the Algorithm 6 of NUTS.
#' @param ... Other parameters that may be useful.
#'
#' @export
NUTS <- function(A, j, iter, lpd_func, grad_lpd_func, M_adapt, M_diag, delta = 0.6, max_treedepth = 10,
                 eps = 1, eps_bar, H, mu, verbose = TRUE, Delta_max, ...){
  kappa <- 0.75
  t0 <- 10
  gamma <- 0.05
  accept <- 0
  dims <- dim(A)

  if(is.null(M_diag)){
    M_diag <- diag(1, dims[1])
  }

  if(iter == 1){
    eps <- find_reasonable_epsilon(A, j, lpd_func, grad_lpd_func, M_diag, eps = eps, verbose = TRUE, ...)
    mu <- log(10*eps)
    H <- 0
    eps_bar <- 1
  }

  aa0 <- rmvnorm(dim = dims[1], mu = 0, sigma = solve_chol(M_diag))

  u <- runif(1, 0, exp(joint_log_density(A, aa0, lpd_func, M_diag, ...)))
  if(is.nan(u)){
    warning("NUTS: sampled slice u is NaN")
    u <- runif(1, 0, 1e5)
  }
  A_minus <- A
  A_plus <- A
  aa_minus <- aa0
  aa_plus <- aa0
  jj <- 0
  n <- 1
  s <- 1
  if(iter > M_adapt){
    eps <- runif(1, 0.9*eps_bar, 1.1*eps_bar)
  }
  while(s == 1){
    # choose direction {-1, 1}
    direction <- sample(c(-1, 1), 1)
    if(direction == -1){
      temp <- build_tree(A=A_minus, j=j, aa=aa_minus, u=u, v=direction, jj=jj, eps=eps, A00=A, aa0=aa0, lpd_func=lpd_func,
                 grad_lpd_func=grad_lpd_func, M_diag=M_diag, Delta_max = Delta_max,...)
      A_minus <- temp$A_minus
      aa_minus <- temp$aa_minus
    } else{
      temp <- build_tree(A=A_plus, j=j, aa=aa_plus, u=u, v=direction, jj=jj, eps=eps, A00=A, aa0=aa0, lpd_func=lpd_func,
                         grad_lpd_func=grad_lpd_func, M_diag=M_diag, Delta_max = Delta_max,...)
      A_plus <- temp$A_plus
      aa_plus <- temp$aa_plus
    }
    if(is.nan(temp$s)) temp$s <- 0
    if(temp$s == 1){
      if(runif(1) < temp$n / n){
        A <- temp$A
        accept <- 1
      }
    }
    n <- n + temp$n
    s <- check_NUTS(temp$s, A_plus, A_minus, aa_plus, aa_minus)
    jj <- jj + 1
    if(jj > max_treedepth){
      warning("NUTS: Reached max tree depth")
      break
    }
  }
  if(iter <= M_adapt){
    H <- (1 - 1/(iter + t0))*H + 1/(iter + t0) * (delta - temp$alpha / temp$n_alpha)
    log_eps <- mu - sqrt(iter)/gamma * H
    eps_bar <- exp(iter^(-kappa) * log_eps + (1 - iter^(-kappa)) * log(eps_bar))
    eps <- exp(log_eps)
  } else{
    eps <- eps_bar
  }

  lpd <- lpd_func(A=A, ...)

  return(list(A = A, lpd = lpd, accept = accept, eps = eps, eps_bar = eps_bar, H = H, mu = mu,
              M_adapt = M_adapt, M_diag = M_diag))
}
