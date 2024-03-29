% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI.R
\name{CI}
\alias{CI}
\title{Calculation of the posterior credible intervals for regression parameters}
\usage{
CI(SIMP.fit, combined = T, Chain.no = 1, levels = 0.95)
}
\arguments{
\item{SIMP.fit}{The MCMC output from \code{SIMP()} function.}

\item{combined}{Logical. Indicate whether posterior samples (after burn-in) from all chains 
should be combined to calculate the posterior credible interval. If there is only one chain,
then there is no difference by indicating this parameter to be TRUE or FALSE.}

\item{Chain.no}{If combined = FALSE, we need to indicate which chain should be used for the 
calculation of posterior credible interval. Chain.no must be an integer between 1 and the number of chains.}

\item{levels}{Significance level, default to be 0.95.}
}
\description{
Calculation of the posterior credible intervals for regression parameters
}
\examples{
\dontrun{
library(SIMP)
library(Renvlp)
data(wheatprotein) # Load Renvlp package only for wheatprotein dataset.
set.seed(1)
X1C = wheatprotein[, 4:5]
X1D = as.matrix(wheatprotein[, 8], ncol = 1)
X2 = wheatprotein[, 6:7]
Y = wheatprotein[, 1:3]
MC_output <- SIMP(X1C = X1C, X1D = X1D, X2 = X2,
                  Y = Y, dx = 1, dy = 1, n.iter = 1e4)
CI(MC_output)
}
}
