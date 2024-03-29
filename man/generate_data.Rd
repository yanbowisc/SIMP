% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{generate_data}
\alias{generate_data}
\title{Data generating function according to SIMP}
\usage{
generate_data(
  muX1C.tru,
  muY.tru,
  beta1C.tru,
  beta1D.tru,
  beta2.tru,
  gamma.tru,
  SigmaCD.tru,
  SigmaYX.tru,
  K,
  mu2,
  SigmaX2,
  n,
  r,
  pc,
  pd,
  p2,
  ...
)
}
\arguments{
\item{muX1C.tru}{a vector of length p_C. The true value of mu1C.}

\item{muY.tru}{a vector of length r. The true value of muY.}

\item{beta1C.tru}{a p_C by r matrix. The true value of beta1C.}

\item{beta1D.tru}{a p_D by r matrix. The true value of beta1D.}

\item{beta2.tru}{a p_2 by r matrix. The true value of beta2.}

\item{gamma.tru}{a p_D by p_C matrix. The true value of gamma.}

\item{SigmaCD.tru}{a p_C by p_C matrix. The true value of SigmaCD.}

\item{SigmaYX.tru}{a r by r matrix. The true value of SigmaYX.}

\item{K}{a positive integer. The K in the discrete uniform{0,1,...,K-1} for the generation of X1D.}

\item{mu2}{a vector of length p_2. The true value of mean of Normal distribution for the generation of X2.}

\item{SigmaX2}{a p_2 by p_2 matrix. The true value of Covariance matrix of Normal distribution for the generation of X2.}

\item{n}{Sample size.}

\item{r}{Dimension of response Y.}

\item{pc}{Dimension of X1C, the continuous part of the predictors of interest.}

\item{pd}{Dimension of X1D, the discrete part of the predictors of interest.}

\item{p2}{Dimension of X2, predictors of not main interest.}

\item{...}{Other parameters needed}
}
\description{
Data generating function according to SIMP
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
 dat <- do.call(generate_data, c(all_pars, 
 list(K = K, mu2 = mu2, SigmaX2 = SigmaX2, n = n, r = r, pc = pc, pd = pd, p2 = p2)))
}
}
