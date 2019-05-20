#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_minn_diffuse <- function(x, ...){

  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(lambda1 = x$lambda1, lambda2 = x$lambda2, lambda3 = x$lambda3,
                                           prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
                                    n_lags = x$n_lags, block_exo = x$block_exo)
  prior_Pi_mean <- matrix(0, n_vars, n_vars*x$n_lags + 1)
  prior_Pi_mean[, 2:(n_vars+1)] <- diag(x$prior_Pi_AR1)

  Y <- x$Y
  freq <- x$freq
  n_fcst <- x$n_fcst
  verbose <- x$verbose
  n_lags <- x$n_lags
  lambda4 <- x$lambda4


  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_Z <- init$init_Z

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  n_q <- sum(freq == "q")
  if (n_q < n_vars) {
    T_b <- max(which(!apply(apply(Y[, freq == "m", drop = FALSE], 2, is.na), 1, any)))
  } else {
    T_b <- nrow(Y)
  }
  if (n_q > 0) {
    Lambda_ <- mfbvar:::build_Lambda(rep("q", n_q), 3)
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
  }

  n_pseudolags <- max(c(n_lags, 3))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags
  d <- matrix(1, nrow = nrow(Y), ncol = 1)

  ################################################################
  ### Preallocation
  # Pi and Sigma store their i-th draws in the third dimension, psi
  # is vectorized so it has its i-th draw stored in the i-th row
  # Pi:    p * pk * n_reps, each [,,i] stores Pi'
  # Sigma: p * p  * n_reps
  # psi:   n_reps * p
  # Z:     T * p * n_reps
  ### If forecasting (h is horizon):
  # Z_fcst: hk * p * n_reps
  # d_fcst_lags: hk * m
  ### If root checking:
  # roots: n_reps vector
  # num_tries: n_reps vector
  ### If smoothing of the state vector:
  # smoothed_Z: T * p * n_reps

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags + 1, n_reps/n_thin))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps/n_thin))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps/n_thin))

  Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }




  ################################################################
  ### MCMC sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the MCMC sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  if (is.null(init_Z)) {
    Z[,, 1] <- mfbvar:::fill_na(Y)
  } else {
    if (all(dim(Z[,, 1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,, 1]), collapse = " x ")))
    }

  }

  ols_results <- mfbvar:::ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = 1)

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- cbind(ols_results$const, ols_results$Pi)
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    if (all(dim(Sigma[,,1]) == dim(init_Sigma))) {
      Sigma[,, 1] <- init_Sigma
    } else {
      stop(paste0("The dimension of init_Sigma is ", paste(dim(init_Sigma), collapse = " x "), ", but should be ", paste(dim(Sigma[,,1]), collapse = " x ")))
    }
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]

  # For the posterior of Pi
  inv_prior_Pi_Omega <- diag(1/c(prior_Pi_Omega))
  Omega_Pi <- matrix(inv_prior_Pi_Omega %*% c(prior_Pi_mean), n_vars*n_lags + 1, n_vars)

  mfbvar:::mcmc_minn_diffuse(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,Lambda_,prior_Pi_Omega,
                             Omega_Pi,Z_1,n_reps,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,
                        n_thin,verbose)


  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = NULL, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, smoothed_Z = NULL, n_determ = 1,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = NULL, prior_psi_mean = NULL, n_reps = n_reps, Lambda_ = Lambda_, freq = freq,
                     init = list(init_Pi = Pi[,, n_reps/n_thin], init_Sigma = Sigma[,, n_reps/n_thin], init_Z = Z[,, n_reps/n_thin]))

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}
