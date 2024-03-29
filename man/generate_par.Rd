% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_par.R
\name{generate_par}
\alias{generate_par}
\title{Parameter generating function in Simulation section}
\usage{
generate_par(r, pc, pd, p2, dx, dy)
}
\arguments{
\item{r}{Dimension of response Y.}

\item{pc}{Dimension of X1C, the continuous part of the predictors of interest.}

\item{pd}{Dimension of X1D, the discrete part of the predictors of interest.}

\item{p2}{Dimension of X2, predictors of not main interest.}

\item{dx}{Partial predictor envelope dimension of SIMP.}

\item{dy}{Partial response envelope dimension of SIMP.}
}
\description{
Parameter generating function in Simulation section
}
\examples{
\dontrun{
r <- 8
pc <- 8
pd <- 2
p2 <- 2
p = pc + pd + p2
K <- 3
mu2 <- c(2, 5)
dx.tru <- 6
dy.tru <- 2
n <- 300
set.seed(2)
if (p2 > 0){
  SigmaX2 <- rinvwish(p2, diag(1, p2), p2)
}else{
  SigmaX2 <- 0
}
all_pars <- generate_par(r, pc, pd, p2, dx.tru, dy.tru)
}
}
