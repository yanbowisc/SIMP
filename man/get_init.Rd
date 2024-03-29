% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_init.R
\name{get_init}
\alias{get_init}
\title{Warm-start estimator for MCMC algorithm of SIMP}
\usage{
get_init(X1C, X1D_ctr, X2_ctr, Y, dx, dy, method.idx = 1)
}
\arguments{
\item{X1C}{Design matrix of the continuous part of the predictors of interest.}

\item{X1D_ctr}{Centered design matrix of the discrete part of the predictors of interest.}

\item{X2_ctr}{Centered design matrix of the predictors of not main interest}

\item{Y}{Response data matrix.}

\item{dx}{Partial predictor envelope dimension. Must be an integer between 0 and p_C.}

\item{dy}{Partial response envelope dimension. Must be an integer between 0 and r.}

\item{method.idx}{Index for the warm start method, default to be 1.}
}
\description{
Warm-start estimator for MCMC algorithm of SIMP
}
