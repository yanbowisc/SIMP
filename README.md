# SIMP
This is a Github repository for `SIMP`, the Bayesian simultaneous partial envelope model. Compared with the simultaneous envelope and the partial response envelope models, this proposed model addresses some limitations of these two advanced envelope models by partitioning the predictors of interest into the continuous and discrete parts. It also provides a more efficient estimator for the regression coefficients, and improves prediction for the responses in the multivariate linear regression model. See details in the manuscript (Shen et al.), which is under review.

## Install `SIMP` R package: 

```R
install.packages("devtools")
devtools::install_github("yanbowisc/SIMP")
```

## Readme description for producing the simulation results of SIMP with r = 3 and n = 300 (The first result of SIMP in Table 3) in Section 7.2.2 of the paper "Bayesian simultaneous partial envelope model with application to an imaging genetics analysis". R codes for all other tables and figures can be similarly obtained.

```R
library(SIMP)
```
`Load the library of SIMP`



