# SIMP
This is a Github repository for `SIMP`, the Bayesian simultaneous partial envelope model. Compared with the simultaneous envelope and the partial response envelope models, this proposed model addresses some limitations of these two advanced envelope models by partitioning the predictors of interest into the continuous and discrete parts. It also provides a more efficient estimator for the regression coefficients, and improves prediction for the responses in the multivariate linear regression model. See details in the manuscript (Shen et al.), which is under review.

## Install `SIMP` R package: 

```R
install.packages("devtools")
devtools::install_github("yanbowisc/SIMP")
```

## Readme description for the demo of using our package

Here, we produce one replicate of the simulation results of SIMP with r = 3 and n = 300 (The first result of SIMP in Table 3) in Section 7.2.2 of the paper "Bayesian simultaneous partial envelope model with application to an imaging genetics analysis". R codes for all other tables and figures can be similarly obtained.

```R
library(SIMP)
# Load the library of SIMP

rep <- 1
# This is the first replication. 

n <- 300
# The sample size is 300

r <- 3
# The dimension of Y is 3

pc <- 8
# The dimension of X1C is 8

pd <- 2
# The dimension of X1D is 2

p2 <- 2
# The dimension of X2 is 2

K <- 2
# Number of levels for generating uniform X1D

mu2 <- c(2, 5)
# means for X2

burnin.prop <- 0.5
# Burn-in proportion for the MCMC algorithm.

dx.tru <- 2
# True envelope dimension dx is set as 2.

dy.tru <- 2
# True envelope dimension dy is set as 2.

method.idx <- 1
# The method of initialization for SIMP is chosen as the first method, our default method.

n.iter <- 2e4
# The number of total steps in our MCMC algorithm is set as 20,000.
```
We load the `SIMP` package and set some parameters in the above codes.

```R
set.seed(1)

if (p2 > 0){

  SigmaX2 <- rinvwish(p2, diag(1, p2), p2)
  
}else{

  SigmaX2 <- 0
  
}
# Generate the covariance parameter for generating X2.

all_pars <- generate_par(r, pc, pd, p2, dx = dx.tru, dy = dy.tru)
# Generate other parameters
```
We generate all parameters that are needed for data generation in the above codes.

```R
set.seed(rep)  
# Set the random seed. It should be adjusted to different values on different replications.

dat <- do.call(generate_data, c(all_pars, list(K = K, mu2 = mu2, SigmaX2 = SigmaX2, n = n, r = r, pc = pc, pd = pd, p2 = p2)))
# Generate the data matrices Y, X1C, X1D and X2.

X1C <- dat$X1C
# Extract the data matrix X1C, with dimensions n * pc.

X1D <- dat$X1D
# Extract the data matrix X1D, with dimensions n * pd.

X2 <- dat$X2
# Extract the data matrix X2, with dimensions n * p2.

X <- dat$X
# Extract the combined predictor matrix X = (X1C, X1D, X2), with dimensions n * (pc + pd + p2).

Y <- dat$Y
# Extract the data matrix Y, with dimensions n * r. 

X1C_bar <- colMeans( X1C )
X1C_ctr <- sweep(X1C, 2, X1C_bar, "-")
# Standardize the data matrix X1C, by subtracting the corresponding mean for each column. The standardized X1C is X1C_ctr.
  
if (pd > 0){
  X1D_bar <- colMeans( X1D )
  X1D_ctr <- sweep(X1D, 2, X1D_bar, "-")
}else{
  X1D_bar <- NULL
  X1D_ctr <- NULL
}
# Standardize the data matrix X1D, by subtracting the corresponding mean for each column. The standardized X1D is X1D_ctr.

X1c <- cbind(X1C_ctr, X1D_ctr)
# Combine the standardized X1C_ctr and X1D_ctr, get the standardized predictor of interest as X1c.
  
if (p2 > 0){
  X2_bar <- colMeans(X2)
  X2_ctr <- sweep(X2, 2, X2_bar, "-")
}else{
  X2_bar <- NULL
  X2_ctr <- NULL
}
# Standardize the data matrix X2, by subtracting the corresponding mean for each column. The standardized X2 is X2_ctr.
  
X_ctr <- cbind(X1c, X2_ctr)
# Combine X1c and X2_ctr to get the standardized full predictor data matrix.

X_std <- scale(X_ctr)
X_sd <- attr(X_std, "scaled:scale")
# X_ctr is further standarized by dividing the corresponding standard deviation for each column, and get X_std.
# The standard deviation for all columns are saved in X_sd.

Y_bar <- colMeans(Y)
Y_ctr <- sweep(Y, 2, Y_bar, "-")
# Standardize the data matrix Y, by subtracting the corresponding mean for each column. The standardized Y is Y_ctr.
```
We generate all data matrices X1C (n * pc), X1D  (n * pd), X2  (n * p2) and Y  (n * r), and standardize them in the above codes. Notice, the random seed 'rep' should be adjusted to different values in different replications (we have 500 replications in total).

```R
bic <- numeric((pc + 1) * (r + 1))
# initialized an empty vector to save BIC-MCMC values for all envelope dimensions

dx_min <- dy_min <- Bura_Cook(X1C_ctr, X1D_ctr, X2_ctr, Y_ctr, sig.level = 0.05)
# Select the minimum envelope dimensions dx_min for dx, and dy_min for dy by the Bura-Cook estimator that is introduced in Section 6.1 in the paper.

for (dx in dx_min:pc){

  for (dy in dy_min:r){
    
    SIMP.fit <- SIMP(X1C, X1D, X2, Y, 
                     dx = dx, 
                     dy = dy,
                     n.iter = n.iter,
                     n.chains = 1,
                     tau = 0.1,
                     init_method = "envlps",
                     Metro_method = "RW",
                     HMC_steps = NA,
                     autotune = TRUE,
                     tune.accpt.prop.lower = 0.3,
                     tune.accpt.prop.upper = 0.4,
                     tune.incr = 0.05,
                     delta = delta,
                     burnin.prop = burnin.prop,
                     tune.burnin.prop = 0.5,
                     tune_nterm = 50,
                     show_progress = TRUE,
                     chains_parallel = FALSE,
                     method.idx = method.idx)
    
    bic[dx*(r + 1) + dy + 1] <- BIC(SIMP.fit, r, pc, pd, p2, dx, dy, n.iter, burnin.prop, samp.size)
    
    
  }
}
# Fit SIMP for each tentative (dx, dy) combination, and save their BIC-MCMC values in the bic vector.

dx_star_SIMP <- ceiling( which( bic == min( bic[ bic != 0 ])) / ( r + 1)) - 1
# save the selected dx as dx_star_SIMP.

dy_star_SIMP <- ifelse(( which( bic == min( bic[ bic != 0 ])) %% ( r + 1) - 1) >= 0, which( bic == min( bic[ bic != 0])) %% ( r + 1 ) - 1, r)
# save the selected dy as dy_star_SIMP.
    
```
We implement the model selection by two steps in the above codes.

(1) We narrow down the range of dimensions to be searched by the Bura-Cook estimator that is introduced in Section 6.1 in the paper.

(2) We select the envelope dimension dx and dy by BIC-MCMC in the narrowed searching range.

```R
SIMP.fit <- SIMP(X1C, X1D, X2, Y, 
                     dx = dx_star_SIMP, 
                     dy = dy_star_SIMP,
                     n.iter = n.iter,
                     n.chains = 1,
                     tau = 0.1,
                     init_method = "envlps",
                     Metro_method = "RW",
                     HMC_steps = NA,
                     autotune = TRUE,
                     tune.accpt.prop.lower = 0.3,
                     tune.accpt.prop.upper = 0.4,
                     tune.incr = 0.05,
                     delta = delta,
                     burnin.prop = burnin.prop,
                     tune.burnin.prop = 0.5,
                     tune_nterm = 50,
                     show_progress = TRUE,
                     chains_parallel = FALSE,
                     method.idx = method.idx)

# Model fitting on SIMP using the selected envelope dimensions dx_star_SIMP and dy_star_SIMP.                   

muY.est <- elmwise_mean_in_list(SIMP.fit$muY[[1]])
mu1C.est <- elmwise_mean_in_list(SIMP.fit$muX1C[[1]])
beta1C.est <- elmwise_mean_in_list(SIMP.fit$beta1C[[1]])
beta1D.est <- elmwise_mean_in_list(SIMP.fit$beta1D[[1]])
beta2.est <- elmwise_mean_in_list(SIMP.fit$beta2[[1]])
# Extract the posterior mean estimators of muY, mu1C, beta1C, beta1D and beta2 from SIMP
```
We implement the model fitting of SIMP with selected envelope dimensions, and calculate the posterior mean estimators of muY, mu1C, beta1C, beta1D and beta2 in the above codes. In our simulation in Section 7.2.2, we repeated the data generation, model selection and model fitting for 500 times, and reported the MSEs of the estimated beta1C and beta1D over 500 replications in Table 3.
