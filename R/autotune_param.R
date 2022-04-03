#' Compute the means of each entry in a matrix over a list
#' @param lst A list of matrices of the same size.
#'
#' @export
elmwise_mean_in_list <- function(lst) {
  n_lst <- length(lst)
  Reduce('+', lst) / n_lst
}

#' Compute the standard deviations of each entry in a matrix over a list
#' @param lst A list of matrices of the same size.
#'
#' @export
elmwise_sd_in_list <- function(lst) {
  n_lst <- length(lst)
  xbar <- elmwise_mean_in_list(lst)
  x2bar <- elmwise_mean_in_list(
    lapply(lst, "^", 2)
  )
  sqrt(abs(x2bar - xbar^2) * n_lst/(n_lst - 1))
}

#' The function to control the acceptance rate of the Metropolis steps by tuning tau.
#' @param draw_param_list The list of all iterates of A or B till the current iteration.
#' @param tune_param The matrix of tau, i.e. the standard deviations of the proposal distribution of A or B.
#' @param accpt_list The list of the acceptance for all iterates of A or B till the current iteration.
#' @param tune_nterm The iteration after which we start to tune tau.
#' @param tune.incr The adjustment proportion for tuning tau.
#' @param tune.accpt.prop.lower The ideal lower bound for the acceptance rate of the Metropolis steps.
#' @param tune.accpt.prop.upper The ideal upper bound for the acceptance rate of the Metropolis steps.
#' @param ... Other parameters that may be useful.
#'
#' @export
autotune_param <- function(
  draw_param_list,
  tune_param,
  accpt_list,
  tune_nterm,
  tune.incr,
  tune.accpt.prop.lower,
  tune.accpt.prop.upper,
  ...
) {

  n_full <- length(draw_param_list)

  ave_accpt_last_nterm <- elmwise_mean_in_list(
    accpt_list[(n_full-tune_nterm+1):n_full]
  )

  idx_increase <- which(ave_accpt_last_nterm > tune.accpt.prop.upper)
  idx_decrease <- which(ave_accpt_last_nterm < tune.accpt.prop.lower)

  if (length(idx_increase) > 0) {
    tune_param[idx_increase] <- tune_param[idx_increase] * (1 + tune.incr)
  }

  if (length(idx_decrease) > 0) {
    tune_param[idx_decrease] <- tune_param[idx_decrease] * (1 - tune.incr)
  }

  tune_param
}
