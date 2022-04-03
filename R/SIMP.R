#' Fit Bayesian simultaneous partial envelope model
#'
#' @import stats parallel methods Renvlp madness
#' @param X1C Design matrix of the continuous part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X1D Design matrix of the discrete part of the predictors of interest.
#' Must have the same number of rows as \code{Y}.
#' @param X2 Design matrix of the nuisance predictors. Must have the same number of rows as \code{Y}.
#' @param Y Response matrix. Must have the same number of rows as \code{X}.
#' @param dx Partial predictor envelope dimension. Must be an integer between 0 and \code{ncol(X1C)}.
#' @param dy Partial response envelope dimension. Must be an integer between 0 and \code{ncol(Y)}.
#' @param n.iter Number of Markov chain iterations to run in each chains. *Includes burn-in*.
#' @param n.chains Number of independent chains to run.
#' @param tau The Metropolis tuning parameter
#' @param tune.accpt.prop.lower Lower bound for the acceptance rate of the metropolis step
#' @param tune.accpt.prop.upper Upper bound for the acceptance rate of the metropolis step
#' @param tune.incr Adjustment magnitude in tuning tau.
#' @param autotune logical. Should the Metropolis tuning parameter
#'  be tuned during burn-in via adaptation?
#' @param init_method Initialization. Available options are "envlps" for the consistent initialization
#' we proposed , "generate" for random generation, or "input"
#' for giving parameters via input param init_params.
#' @param init_params Input parameter values if the initial method is "input".
#' @param Metro_method method for metropolis step for sampling A or B. Available options are "RW" for
#' Random walk metropolis, "HMC" for Hamiltonian monte carlo, and "NUTS" for No-U-Turn sampler.
#' "RW" is recommended to use only for accuracy and especially efficiency through our testing.
#' @param HMC_steps Number of steps in Hamiltonian monte carlo in each iteration.
#' @param burnin.prop Proportion for burn-in period among all iterations.
#' @param tune.burnin.prop Proportion for the Metropolis tuning parameter be tuned
#' among burn-in period.
#' @param tune_nterm After which iteration the Metropolis tuning parameter should be tuned.
#' @param show_progress Logical. Indicate whether the progress bar show or off.
#' @param chains_parallel Logical. Indicate whether we should use parallel computing for each chain.
#' @param cores Number of cores used in parallel computing.
#' @param method.idx Choice of the method of initialization. There are 6 methods in total. The default one is
#' the first one (consistent estimator in the manuscript). Always keep to be the default value
#' if no special reasons.
#' @param random.seed Whether we should fix random seed at the start of running for each chain.
#' @param ... all other useful inputs.
#' @references Yanbo Shen, Yeonhee Park, Saptarshi Chakraborty, Chunming Zhang(202X)
#' @examples
#' \dontrun{
#' library(SIMP)
#' library(Renvlp)
#' data(wheatprotein) # Load Renvlp package only for wheatprotein dataset.
#' set.seed(1)
#' X1C = wheatprotein[, 4:5]
#' X1D = as.matrix(wheatprotein[, 8], ncol = 1)
#' X2 = wheatprotein[, 6:7]
#' Y = wheatprotein[, 1:3]
#' MC_output <- SIMP(X1C = X1C, X1D = X1D, X2 = X2,
#'                   Y = Y, dx = 1, dy = 1, n.iter = 1e4)
#' }
#' @export
SIMP <- function(X1C, X1D, X2, Y, dx, dy,
                 n.iter = 2e4,
                 n.chains = 1,
                 tau = 0.1,
                 init_method = "envlps",
                 init_params = NULL,
                 Metro_method = "RW",
                 HMC_steps = 10,
                 autotune = TRUE,
                 tune.accpt.prop.lower = 0.3,
                 tune.accpt.prop.upper = 0.4,
                 tune.incr = 0.05,
                 burnin.prop = 0.5,
                 tune.burnin.prop = 0.5,
                 tune_nterm = 50,
                 show_progress = TRUE,
                 chains_parallel = FALSE,
                 cores = 1,
                 method.idx = 1,
                 random.seed = T,
                 ...)
{
  # Data pre-processing
  dim_Y <- dim(Y)
  n <- dim_Y[1]
  r <- dim_Y[2]
  pc <- ncol(X1C)
  if (is.null(X1D)){
    pd <- 0
  }else{
    pd <- ncol(X1D)
  }

  if (is.null(X2)){
    p2 <- 0
  }else{
    p2 <- ncol(X2)
  }

  #-------------------------------------------------------------------------------------
  # Function parameters checking
  #-------------------------------------------------------------------------------------
  if (dx < 0 || dx > pc || dx != as.integer(dx))
    stop(paste0("\'dx\' must be an integer between 0 and ", pc, " (pc)") )
  if (dy < 0 || dy > r || dy != as.integer(dy))
    stop(paste0("\'dy\' must be an integer between 0 and ", r, " (r)") )

  Y <- Y.orig <- data.matrix(Y)
  Y_bar <- colMeans(Y)

  X1C <- X1C.orig <- data.matrix(X1C)
  X1C_bar <- colMeans(X1C)
  Z_bar <- c( X1C_bar, Y_bar)

  if (pd > 0){
    X1D <- X1D.orig <- data.matrix(X1D)
    X1D_bar <- colMeans(X1D)
    X1D_ctr <- sweep(X1D, 2, X1D_bar, "-")
    X1D_ctr.t_X1D_ctr <- crossprod(X1D_ctr)
  }else{
    X1D <- X1D.orig <- X1D_bar <- X1D_ctr <- X1D_ctr.t_X1D_ctr <- NULL
  }

  if (p2 > 0){
    X2 <- X2.orig <- data.matrix(X2)
    X2_bar <- colMeans(X2)
    X2_ctr <- sweep(X2, 2, X2_bar, "-")
    X2_ctr.t_X2_ctr <- crossprod(X2_ctr)
  }else{
    X2 <- X2.orig <- X2_bar <- X2_ctr <- X2_ctr.t_X2_ctr <- NULL
  }


  # Metropolis parameter setup

  n.burnin <- ceiling(n.iter * burnin.prop)
  burnin_iter <- seq_len(n.burnin)
  n.tune <- ceiling(n.burnin * tune.burnin.prop)
  tune_iter <- seq_len(n.tune)

  # Set Hyper-parameters
  #hyper_params_priors <- list(e = matrix(rep(0, 6), 2, 3))
  hp <- get_hyperparams_init(pc, pd, p2, r, dx, dy)
  # parameters initialization, for simplicity use true values first.
  if(init_method == "envlps"){
    params_init <- get_init(X1C, X1D_ctr, X2_ctr, Y, dx, dy, method.idx = method.idx)
  }else if(init_method == "generate"){
    params_init <- generate_par(r, pc, pd, p2, dx, dy)
    names(params_init) <- c("muX1C", "muY", "beta1C", "beta1D", "beta2", "gamma", "Omega", "Omega0",
                            "Phi", "Phi0", "A", "L", "L0","B", "R", "R0","etaC", "etaD", "SigmaCD", "SigmaYX")
  }else if(init_method == "input"){
    params_init <- init_params
  }

  Yctr <- params_init$Yc
  Yctr.t_Yctr <- crossprod(Yctr)
  X1C_ctr <- params_init$X1Cc
  X1C_ctr.t_X1c_ctr <- crossprod(X1C_ctr)
  X.order <- params_init$X.order
  Y.order <- params_init$Y.order
  X1C <- X1C[, X.order]
  Y <- Y[, Y.order]
  init_idx <- params_init$init_idx

  runMC <- function(params_init = NULL, chain_no){
    #-----------------------------runMC starts-----------------------------------------------------
    # takeout initial params and prepare for MCMC iterations
    A_exists <- (dx > 0) & (dx < pc)
    B_exists <- (dy > 0) & (dy < r)
    one_n <- rep(1, n)
    muX1C <- params_init$muX1C
    muY <- params_init$muY
    beta1C <- params_init$beta1C
    beta1D <- params_init$beta1D
    beta2 <- params_init$beta2
    gamma <- params_init$gamma
    etaC <- params_init$etaC
    etaD <- params_init$etaD
    if (A_exists){
      A <- params_init$A
      L_L0 <- find_gammas_from_A(A)
      L <- L_L0$gamma
      L0 <- L_L0$gamma0
    }else if(dx == 0){
      A <- 0
      L <- 0
      L0 <- diag(1, pc)
    }else if(dx == pc){
      A <- 0
      L0 <- 0
      L <- diag(1, pc)
    }

    if (B_exists){
      B <- params_init$B
      R_R0 <- find_gammas_from_A(B)
      R <- R_R0$gamma
      R0 <- R_R0$gamma0
    }else if(dy == 0){
      B <- 0
      R <- 0
      R0 <- diag(1, r)
    }else if(dy == r){
      B <- 0
      R0 <- 0
      R <- diag(1, r)
    }


    Omega <- params_init$Omega
    if (dx > 0)  {
      Omega.inv <- solve_chol(Omega)
      Omega.half <- sqrtmat(Omega)
    } else {
      Omega.inv <- 0
      Omega.half <- 0
    }

    Omega0 <- params_init$Omega0
    if (dx < pc)  {
      Omega0.inv <- solve_chol(Omega0)
    } else {
      Omega0.inv <- 0
    }

    Phi <- params_init$Phi
    if (dy > 0)  {
      Phi.inv <- solve_chol(Phi)
    } else {
      Phi.inv <- 0
      Phi.half.inv <- 0
    }

    Phi0 <- params_init$Phi0
    if (dy < r)  {
      Phi0.inv <- solve_chol(Phi0)
    } else {
      Phi0.inv <- 0
    }

    if (A_exists) {
      SigmaCD <- L %*% tcrossprod(Omega, L) +
        L0 %*% tcrossprod(Omega0, L0)
    } else if (dx == 0) {
      SigmaCD <- Omega0
    } else if (dx == pc) {
      SigmaCD <- Omega
    }
    SigmaCD.inv <- solve_chol(SigmaCD)

    if (B_exists) {
      SigmaYX <- R %*% tcrossprod(Phi, R) +
        R0 %*% tcrossprod(Phi0, R0)
    } else if (dy == 0) {
      SigmaYX <- Phi0
    } else if (dy == r) {
      SigmaYX <- Phi
    }
    SigmaYX.inv <- solve_chol(SigmaYX)
    if (pd > 0){
      Lambda.inv <- solve_chol(hp$Lambda)
      Lambda.half <- sqrtmat(hp$Lambda)
      Q.half <- sqrtmat(hp$Q)
      Q.inv <- solve_chol(hp$Q)
    }else{
      Lambda.inv <- Lambda.half <- Q.half <-  Q.inv <- 0
    }

    if (p2 > 0){
      M.half <- sqrtmat(hp$M)
      M.inv <- solve_chol(hp$M)
      beta2.shifted <- beta2 - M.inv%*%hp$e
    }else{
      M.half <- 0
      M.inv <- 0
      beta2.shifted <- 0
    }


    if (dx > 0){
      E.half <- sqrtmat(hp$E)
      E.inv <- solve_chol(hp$E)
    }else{
      E.half <- 0
      E.inv <- 0
    }


    if (dy * pd > 0){
      X1D_ctr_etaD.t <- tcrossprod(X1D_ctr, etaD)
      etaD_shifted <- etaD - hp$f%*%Q.inv
    }else{
      X1D_ctr_etaD.t <- etaD_shifted <- 0
    }

    if (A_exists){
      tau.A <- matrix(tau, pc - dx, dx)
      eps.A <- eps_bar.A <- H.A <- mu.A <- matrix(1, pc - dx, dx)
    }else{
      tau.A <- 0
    }

    if (B_exists){
      tau.B <- matrix(tau, r - dy, dy)
      eps.B <- eps_bar.B <- H.B <- mu.B <- matrix(1, r - dy, dy)
    }else{
      tau.B <- 0
    }

    if (p2 == 0){
      X2_ctr_beta2 <- 0
    }

    if (pd == 0){
      gamma_shifted <- 0
    }

    if(dx == 0){
      X1C_ctr_L <- 0
    }

    if (dx*dy == 0){
      X1C_ctr_L_etaC.t <- etaC_shifted <- 0
    }

    DeltaZ <- diag(0, pc + r) # just initializes, will be replaced at every iteration

    beta2.list <- gamma.list <- muZ.list <- muX1C.list <- muY.list <-
      etaC.list <- etaD.list <- Omega.list <-  Omega0.list <-
      Phi.list <-  Phi0.list <-
      A.list <- accpt.A.list <- L.list <- L0.list <-
      B.list <- accpt.B.list <- R.list <- R0.list <-
      beta1C.list <- beta1D.list <- resi.list <-
      SigmaCD.list <- SigmaYX.list <- vector("list", n.iter)

    lpd.full.all <- lpd.A.all <- lpd.B.all <- llik.all <-  accpt.A.ave <- rep(0, n.iter)

    if (show_progress){
      pb <- tryCatch(tkProgressBar(title = paste("Chain", chain_no),
                                   label = "Progress: 0%",
                                   min = 0, max = n.iter),
                     error = function(e) e,
                     warning = function(w) w)

      if (is(pb, "error") | is(pb, "warning")) {
        show_progress <- FALSE
        cat(paste("\'tcltk\' could not be loaded. \'show_progress\' set to FALSE."))
      }
    }

    if (random.seed) set.seed(chain_no)


    for(iter in 1:n.iter){

      # generate beta2, if existed (p2>0)
      if(p2 >0){

        M.tilde <- X2_ctr.t_X2_ctr + hp$M
        M.tilde.inv <- solve_chol(M.tilde)
        if (pd > 0){
          e.tilde <- crossprod(X2_ctr, Yctr - X1C_ctr%*%beta1C - X1D_ctr%*%beta1D) + hp$e
        }else{
          e.tilde <- crossprod(X2_ctr, Yctr - X1C_ctr%*%beta1C) + hp$e
        }

        beta2.list[[iter]] <- beta2 <- rMatrixNormal(
          M.tilde.inv%*%e.tilde,
          M.tilde.inv,
          SigmaYX)
        beta2.shifted <- beta2 - M.inv%*%hp$e
        X2_ctr_beta2 <- X2_ctr%*%beta2
      }

      # generate gamma, if existed (pc, pd >0)
      if (pd > 0){
        Lambda.tilde <- X1D_ctr.t_X1D_ctr + hp$Lambda
        Lambda.tilde.inv <- solve_chol(Lambda.tilde)
        g.tilde <- crossprod(X1D_ctr, X1C_ctr) + hp$g
        gamma.list[[iter]] <- gamma <- rMatrixNormal(
          Lambda.tilde.inv%*%g.tilde,
          Lambda.tilde.inv,
          SigmaCD)
        gamma_shifted <- gamma - Lambda.inv%*%hp$g
      }


      # generate muX1C and muY, if existed (pc > 0 (r>0 is assumed to be always true))

      DeltaZ[1:pc, 1:pc] <- SigmaCD
      if(dx * dy > 0){
        DeltaZ[(pc + 1):(pc + r), (pc + 1):(pc + r)] <- SigmaYX + tcrossprod(R%*%etaC%*%Omega.half)
        DeltaZ[1:pc, (pc + 1):(pc + r)] <- tcrossprod(L%*%Omega, R%*%etaC)
        DeltaZ[(pc + 1):(pc + r), 1:pc] <- t(DeltaZ[1:pc, (pc + 1):(pc + r)])
      }else{
        DeltaZ[(pc + 1):(pc + r), (pc + 1):(pc + r)] <- SigmaYX
        DeltaZ[(pc + 1):(pc + r), 1:pc] <- DeltaZ[1:pc, (pc + 1):(pc + r)] <- 0
      }



      muZ.list[[iter]] <- muZ <- rmvnorm(dim = pc + r, mu = Z_bar, sigma = DeltaZ/n)
      muX1C.list[[iter]] <- muX1C <- muZ[1:pc]
      muY.list[[iter]] <- muY <- muZ[(pc+1):(pc+r)]

      X1C_ctr <- sweep(X1C, 2, muX1C, "-")
      Yctr <- sweep(Y, 2, muY, "-")
      if (p2 > 0){
        Yctr_shifted <- Yctr - X2_ctr%*%beta2
      }else{
        Yctr_shifted <- Yctr
      }

      if (pd > 0){
        X1C_ctr_shifted <- X1C_ctr - X1D_ctr%*%gamma
      }else{
        X1C_ctr_shifted <- X1C_ctr
      }

      if (dx > 0){
        X1C_ctr_L <- X1C_ctr%*%L
      }



      # generate etaC, if existed (dx, dy >0)
      if (dx * dy >0){
        E.tilde <- crossprod(X1C_ctr_L) + hp$E
        E.tilde.inv <- solve_chol(E.tilde)
        h.tilde <- crossprod(Yctr_shifted%*%R - X1D_ctr_etaD.t, X1C_ctr_L) + hp$h
        etaC.list[[iter]] <- etaC <- rMatrixNormal(
          h.tilde%*%E.tilde.inv,
          Phi,
          E.tilde.inv
        )
        X1C_ctr_L_etaC.t <- tcrossprod(X1C_ctr_L, etaC)
        etaC_shifted <- etaC - hp$h%*%E.inv
      }


      # generate etaD, if existed (pd, dy > 0)
      if ( pd * dy > 0){
        Q.tilde <- X1D_ctr.t_X1D_ctr + hp$Q
        Q.tilde.inv <- solve_chol(Q.tilde)
        f.tilde <- crossprod( Yctr_shifted%*%R - X1C_ctr_L_etaC.t, X1D_ctr) + hp$f
        etaD.list[[iter]] <- etaD <- rMatrixNormal(
          f.tilde%*%Q.tilde.inv,
          Phi,
          Q.tilde.inv
        )
        X1D_ctr_etaD.t <- tcrossprod(X1D_ctr, etaD)
        etaD_shifted <- etaD - hp$f%*%Q.inv
      }

      # generate Omega, if existed (dx > 0)
      if (dx > 0){
        if (pd > 0){
          PsiX.tilde <- crossprod(X1C_ctr_shifted%*%L) + crossprod(Lambda.half%*%gamma_shifted%*%L) + hp$PsiX
        }else{
          PsiX.tilde <- crossprod(X1C_ctr_shifted%*%L) + hp$PsiX
        }

        wX.tilde <- n + pd + hp$wX
        Omega.list[[iter]] <- Omega <-
          rinvwish(dim = dx, Phi = PsiX.tilde, nu =  wX.tilde)
        Omega.inv <- solve_chol(Omega)
        Omega.half <- sqrtmat(Omega)
      }

      # generate Omega0, if existed (dx < pc)
      if (dx < pc){
        if (pd > 0){
          PsiX0.tilde <- crossprod(X1C_ctr_shifted%*%L0) + crossprod(Lambda.half%*%gamma_shifted%*%L0) + hp$PsiX0
        }else{
          PsiX0.tilde <- crossprod(X1C_ctr_shifted%*%L0) + hp$PsiX0
        }

        wX0.tilde <- n + pd + hp$wX0
        Omega0.list[[iter]] <- Omega0 <-
          rinvwish(dim = pc - dx, Phi = PsiX0.tilde, nu =  wX0.tilde)
        Omega0.inv <- solve_chol(Omega0)
      }


      # generate Phi, if existed (dy > 0)
      if (dy > 0){
        if (p2 * dx > 0){
          if (pd > 0){
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1C_ctr_L_etaC.t - X1D_ctr_etaD.t) +
              crossprod( M.half %*% beta2.shifted %*% R) +
              tcrossprod(etaC_shifted%*%E.half) +
              tcrossprod(etaD_shifted%*%Q.half) + hp$PsiY
          }else{
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1C_ctr_L_etaC.t) +
              crossprod( M.half %*% beta2.shifted %*% R) +
              tcrossprod(etaC_shifted%*%E.half) + hp$PsiY
          }

        }else if ((p2 == 0)&(dx > 0)){
          if (pd > 0){
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1C_ctr_L_etaC.t - X1D_ctr_etaD.t) +
              tcrossprod(etaC_shifted%*%E.half) +
              tcrossprod(etaD_shifted%*%Q.half) + hp$PsiY
          }else{
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1C_ctr_L_etaC.t) +
              tcrossprod(etaC_shifted%*%E.half) + hp$PsiY
          }

        }else if ((p2 > 0)&(dx == 0)){
          if (pd > 0){
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1D_ctr_etaD.t) +
              crossprod( M.half %*% beta2.shifted %*% R) +
              tcrossprod(etaD_shifted%*%Q.half) + hp$PsiY
          }else{
            PsiY.tilde <- crossprod(Yctr_shifted%*%R) +
              crossprod( M.half %*% beta2.shifted %*% R) + hp$PsiY
          }
        }else if ((p2 == 0)&(dx == 0)){
          if (pd > 0){
            PsiY.tilde <- crossprod( Yctr_shifted%*%R - X1D_ctr_etaD.t) +
              tcrossprod(etaD_shifted%*%Q.half) + hp$PsiY
          }else{
            PsiY.tilde <- crossprod(Yctr_shifted%*%R) + hp$PsiY
          }
        }


        wY.tilde <- n + dx + pd + p2 + hp$wY
        Phi.list[[iter]] <- Phi <-
          rinvwish(dim = dy, Phi = PsiY.tilde, nu =  wY.tilde)
        Phi.inv <- solve_chol(Phi)
        Phi.half.inv <- sqrtmatinv(Phi)
      }

      # generate Phi0, if existed (dy < r)
      if (dy < r){
        if (p2 > 0){
          PsiY0.tilde <- crossprod(Yctr_shifted%*%R0) + crossprod( M.half %*% beta2.shifted %*% R0) + hp$PsiY0
        }else{
          PsiY0.tilde <- crossprod(Yctr_shifted%*%R0) + hp$PsiY0
        }

        wY0.tilde <- n + p2 + hp$wY0
        Phi0.list[[iter]] <- Phi0 <-
          rinvwish(dim = r - dy, Phi = PsiY0.tilde, nu =  wY0.tilde)
        Phi0.inv <- solve_chol(Phi0)
      }

      # update SigmaYX.inv before updating A, as SigmaYX.inv is used in updating A.
      if (B_exists) {
        SigmaYX <- R %*% tcrossprod(Phi, R) + R0 %*% tcrossprod(Phi0, R0)
      } else if (dy == 0) {
        SigmaYX <- Phi0
      } else if (dy == r) {
        SigmaYX <- Phi
      }
      SigmaYX.inv <- solve_chol(SigmaYX)

      # generate A, if existed (0 < dx < pc)
      if (A_exists){
        rmh.A <- rmh_colwise_new(A_start = A,
                                 lpd_func = lpd_A,
                                 method = Metro_method,
                                 tau = tau.A,
                                 eps = eps.A,
                                 eps_bar = eps_bar.A,
                                 H = H.A,
                                 mu = mu.A,
                                 LL = HMC_steps,
                                 ColSigma.A.used = NULL,
                                 alpha = 1,
                                 grad_lpd_func,
                                 M_adapt = n.tune,
                                 M_diag = NULL,
                                 delta = NULL,
                                 iter,
                                 samp.size = n,
                                 X1C_ctr,
                                 X1C_ctr_shifted,
                                 Yctr = Yctr,
                                 etaC,
                                 R,
                                 dy,
                                 pd,
                                 X1D_ctr_etaD.t,
                                 X2_ctr_beta2,
                                 Omega.inv,
                                 Omega0.inv,
                                 SigmaYX.inv,
                                 gamma_shifted,
                                 Lambda.half,
                                 A0 = hp$A0,
                                 KA.half.inv = hp$KA.half.inv,
                                 SigmaA.half.inv = hp$SigmaA.half.inv)
        A.list[[iter]] <- A <- rmh.A$A
        accpt.A.list[[iter]] <- tcrossprod(rep(1, pc-dx), rmh.A$accpt)
        lpd.A <- rmh.A$lpd
        lpd.A.all[iter] <- c(lpd.A)
        L_L0 <- attr(lpd.A, "L_L0")
        L.list[[iter]] <- L <- L_L0$gamma
        L0.list[[iter]] <- L0 <- L_L0$gamma0
      } else if (dx == pc){
        L.list[[iter]] <- L <- diag(1, pc)
      } else if (dx == 0){
        L0.list[[iter]] <- L0 <- diag(1, pc)
      }
      if (dx > 0){
        X1C_ctr_L <- X1C_ctr%*%L
        if (dy > 0){
          X1C_ctr_L_etaC.t <- tcrossprod(X1C_ctr_L, etaC)
        }
      }

      # generate B, if existed (0 < dy < r)
      if (B_exists){
        rmh.B <- rmh_colwise_new(A_start = B,
                                 lpd_func = lpd_B,
                                 method = Metro_method,
                                 tau = tau.B,
                                 eps = eps.B,
                                 eps_bar = eps_bar.B,
                                 H = H.B,
                                 mu = mu.B,
                                 LL = HMC_steps,
                                 ColSigma.A.used = NULL,
                                 alpha = 1,
                                 grad_lpd_func,
                                 M_adapt = n.tune,
                                 M_diag = NULL,
                                 delta = NULL,
                                 iter,
                                 samp.size = n,
                                 X1C_ctr,
                                 Yctr = Yctr,
                                 L,
                                 dx,
                                 pd,
                                 etaC,
                                 X1D_ctr_etaD.t,
                                 X2_ctr_beta2,
                                 Phi.inv,
                                 Phi0.inv,
                                 beta2.shifted,
                                 M.half,
                                 B0 = hp$B0,
                                 KB.half.inv = hp$KB.half.inv,
                                 SigmaB.half.inv = hp$SigmaB.half.inv)
        B.list[[iter]] <- B <- rmh.B$A
        accpt.B.list[[iter]] <- tcrossprod(rep(1, r-dy), rmh.B$accpt)
        lpd.B <- rmh.B$lpd
        lpd.B.all[iter] <- c(lpd.B)
        R_R0 <- attr(lpd.B, "R_R0")
        R.list[[iter]] <- R <- R_R0$gamma
        R0.list[[iter]] <- R0 <- R_R0$gamma0
      } else if (dy == r){
        R.list[[iter]] <- R <- diag(1, r)
      } else if (dy == 0){
        R0.list[[iter]] <- R0 <- diag(1, r)
      }

      # update beta1C, if existed (pc > 0)
      if (dx * dy > 0) {
        beta1C.list[[iter]] <- beta1C <- tcrossprod(L, R%*%etaC)
      }else{
        beta1C.list[[iter]] <- beta1C <- matrix(0, nrow = pc, ncol = r)
      }

      # update beta1D, if existed (pd > 0)
      if (pd * dy > 0){
        beta1D.list[[iter]] <- beta1D <- t(R%*%etaD)
      }else if ((pd > 0)&(dy == 0)){
        beta1D.list[[iter]] <- beta1D <- matrix(0, nrow = pd, ncol = r)
      }else{
        beta1D.list[[iter]] <- beta1D <- 0
      }

      # update SigmaCD
      if (A_exists) {
        SigmaCD.list[[iter]] <- SigmaCD <- L %*% tcrossprod(Omega, L) +
          L0 %*% tcrossprod(Omega0, L0)
      } else if (dx == 0) {
        SigmaCD.list[[iter]] <- SigmaCD <- Omega0
      } else if (dx == pc) {
        SigmaCD.list[[iter]] <- SigmaCD <- Omega
      }
      SigmaCD.inv <- solve_chol(SigmaCD)

      # update SigmaYX
      if (B_exists) {
        SigmaYX.list[[iter]] <- SigmaYX <- R %*% tcrossprod(Phi, R) +
          R0 %*% tcrossprod(Phi0, R0)
      } else if (dy == 0) {
        SigmaYX.list[[iter]] <- SigmaYX <- Phi0
      } else if (dy == r) {
        SigmaYX.list[[iter]] <- SigmaYX <- Phi
      }
      SigmaYX.inv <- solve_chol(SigmaYX)

      if (pd * p2 > 0){
        resi.list[[iter]] <- resi <- Yctr - X1C_ctr%*%beta1C - X1D_ctr%*%beta1D - X2_ctr%*%beta2
      }else if ((pd == 0)&(p2 > 0)){
        resi.list[[iter]] <- resi <- Yctr - X1C_ctr%*%beta1C - X2_ctr%*%beta2
      }else if ((pd > 0)&(p2 == 0)){
        resi.list[[iter]] <- resi <- Yctr - X1C_ctr%*%beta1C - X1D_ctr%*%beta1D
      }else{
        resi.list[[iter]] <- resi <- Yctr - X1C_ctr%*%beta1C
      }


      # update log-likelihood
      log_det_Omega <- ifelse(dx > 0, log(det(Omega)), 0)
      log_det_Omega0 <- ifelse(dx < pc, log(det(Omega0)), 0)
      log_det_Phi <- ifelse(dy > 0, log(det(Phi)), 0)
      log_det_Phi0 <- ifelse(dy < r, log(det(Phi0)), 0)

      llik.all[iter] <- llik <- -0.5*(n*(log_det_Omega + log_det_Omega0) +
                                        sum(X1C_ctr_shifted %*% SigmaCD.inv * X1C_ctr_shifted) +
                                        n*(log_det_Phi + log_det_Phi0) +
                                        sum(resi %*% SigmaYX.inv * resi))

      lpd.full.all[iter] <- lpd.full <- llik - 0.5*((pd + hp$wX + dx + 1)*log_det_Omega +
                                                      (pd + hp$wX0 + pc - dx + 1)*log_det_Omega0 +
                                                      (p2 + hp$wY + dx + dy + pd + 1)*log_det_Phi +
                                                      (p2 + hp$wY0 + r - dy + 1)*log_det_Phi0 +
                                                      sum(hp$PsiX*Omega.inv) +
                                                      sum(hp$PsiX0*Omega0.inv) +
                                                      sum(hp$PsiY*Phi.inv) +
                                                      sum(hp$PsiY0*Phi0.inv))
      if (dx * dy > 0){
        lpd.full.all[iter] <- lpd.full - 0.5 * sum(hp$E*crossprod(Phi.half.inv%*%etaC_shifted))
      }

      if (A_exists){
        lpd.full.all[iter] <- lpd.full <- lpd.full - 0.5 * sum((hp$KA.half.inv %*% (A-hp$A0) %*% hp$SigmaA.half.inv)^2)
      }

      if (B_exists){
        lpd.full.all[iter] <- lpd.full <- lpd.full - 0.5 * sum((hp$KB.half.inv %*% (B-hp$B0) %*% hp$SigmaB.half.inv)^2)
      }

      if (pd > 0){
        lpd.full.all[iter] <- lpd.full <- lpd.full - 0.5 * sum(SigmaCD.inv*crossprod(Lambda.half%*%gamma_shifted))
        if (dy > 0){
          lpd.full.all[iter] <- lpd.full <- lpd.full - 0.5 * sum(hp$Q*crossprod(Phi.half.inv%*%etaD_shifted))
        }
      }

      if (p2 > 0){
        lpd.full.all[iter] <- lpd.full <- lpd.full - 0.5 * sum(SigmaYX.inv*crossprod(M.half%*%beta2.shifted))
      }


      if (autotune & iter >= tune_nterm &
          iter %% 5 == 0 & iter <= n.tune) {

        if (A_exists){
          tau.A <- autotune_param(
            draw_param_list = A.list[1:iter],
            tune_param = tau.A,
            accpt_list = accpt.A.list[1:iter],
            tune_nterm = tune_nterm,
            tune.incr = tune.incr,
            tune.accpt.prop.lower = tune.accpt.prop.lower,
            tune.accpt.prop.upper = tune.accpt.prop.upper
          )
        }

        if (B_exists){
          tau.B <- autotune_param(
            draw_param_list = B.list[1:iter],
            tune_param = tau.B,
            accpt_list = accpt.B.list[1:iter],
            tune_nterm = tune_nterm,
            tune.incr = tune.incr,
            tune.accpt.prop.lower = tune.accpt.prop.lower,
            tune.accpt.prop.upper = tune.accpt.prop.upper
          )
        }

      }

      if (show_progress)
        utils::setTxtProgressBar(pb, iter)
    }

    MC <- list(muX1C = muX1C.list[-burnin_iter],
               muY = muY.list[-burnin_iter],
               beta1C = beta1C.list[-burnin_iter],
               beta1D = beta1D.list[-burnin_iter],
               beta2 = beta2.list[-burnin_iter],
               gamma = gamma.list[-burnin_iter],
               Omega = Omega.list[-burnin_iter],
               Omega0 = Omega0.list[-burnin_iter],
               Phi = Phi.list[-burnin_iter],
               Phi0 = Phi0.list[-burnin_iter],
               A = A.list[-burnin_iter],
               L = L.list[-burnin_iter],
               L0 = L0.list[-burnin_iter],
               B = B.list[-burnin_iter],
               R = R.list[-burnin_iter],
               R0 = R0.list[-burnin_iter],
               etaC = etaC.list[-burnin_iter],
               etaD = etaD.list[-burnin_iter],
               SigmaCD = SigmaCD.list[-burnin_iter],
               SigmaYX = SigmaYX.list[-burnin_iter],
               accpt.A.list = accpt.A.list[-burnin_iter],
               accpt.B.list = accpt.B.list[-burnin_iter],
               llik = llik.all[-burnin_iter],
               lpd.full.all = lpd.full.all[-burnin_iter],
               lpd.A.all = lpd.A.all[-burnin_iter],
               lpd.B.all = lpd.B.all[-burnin_iter],
               params_init = params_init,
               tau.A = tau.A,
               tau.B = tau.B,
               X.order = X.order,
               Y.order = Y.order,
               init_idx = init_idx,
               prior_param = hp)

    if (show_progress){
      close(pb)
    }
    MC

  }


  lapply_ <- function(...) {
    if (n.chains == 1 | !chains_parallel) {
      lapply(...)
    } else {
      mclapply(..., mc.cores = cores)
    }
  }

  all_MCs <- lapply_(1:n.chains,
                     function(j) runMC(params_init, chain_no = j))


  vars <- setdiff(names(all_MCs[[1]]), "prior_param")

  out_vars <- lapply(vars,
                     function(x)
                     {
                       out <- lapply(all_MCs, "[[", x)
                       names(out) <- paste0("Chain_", 1:n.chains)
                       out
                     }
  )
  names(out_vars) <- vars

  out_res <- c(list(n.iter = n.iter,
                    X1C = X1C.orig,
                    X1D = X1D.orig,
                    X2 = X2.orig,
                    Y = Y.orig,
                    X.order = X.order,
                    Y.order = Y.order,
                    dx = dx,
                    dy = dy,
                    burnin = n.burnin,
                    tune = n.tune,
                    n.chains = n.chains,
                    prior_param = all_MCs[[1]]$prior_param),
                    out_vars)

  class(out_res) <- c("Benvlp", "SIMP")

  out_res


}
