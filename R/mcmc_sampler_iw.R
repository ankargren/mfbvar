mfbvar_steadystate_iw <- function(x, ssng, ...) {

  envir <- environment()

  if (ssng) {
    required_params <- c("Y", "d", "prior_psi_mean",
                         "n_lags", "n_burnin", "n_reps")
    prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                      "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                      "lambda1", "lambda3", "prior_Pi_AR1", "prior_ng", "s")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                          "n_T", "n_T_", "n_thin", "n_determ", "Z_1")
    params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "Sigma")
  } else {
    required_params <- c("Y", "d", "prior_psi_mean", "prior_psi_Omega",
                         "n_lags", "n_burnin", "n_reps")
    prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                      "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                      "lambda1", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                          "n_T", "n_T_", "n_thin", "n_determ", "Z_1")
    params <- c("Z", "psi", "Pi", "Sigma")
  }


  # Check inputs
  mfbvar:::check_required_params(x, required_params)

  # Assign variables from prior object
  mfbvar:::list_to_variables(x, envir, prior_params)

  # Retrieve some additional variables
  init_vars <- mfbvar:::variable_initialization(Y = Y, freq = freq, freqs = freqs,
                                                n_lags = n_lags, Lambda_ = Lambda_,
                                                n_thin = n_thin)
  mfbvar:::list_to_variables(init_vars, envir, retrieved_params)

  # Prior
  priors <- mfbvar:::prior_Pi_Sigma(lambda1 = lambda1,
                                           lambda2 = lambda3,
                                           prior_Pi_AR1 = prior_Pi_AR1,
                                           Y = Y,
                                           n_lags = n_lags,
                                           prior_nu = n_vars + 2)
  mfbvar:::list_to_variables(priors, envir, "prior_Pi_mean", "prior_Pi_Omega",
                             "prior_S", "inv_prior_Pi_Omega", "Omega_Pi")

  # Initialize ss objects
  init_ss <- mfbvar:::ss_initialization(d = d, d_fcst = d_fcst,
                                        n_T = n_T, n_lags = n_lags, n_fcst = n_fcst)
  mfbvar:::list_to_variables(init_ss, envir, "d_fcst_lags", "D_mat", "dt", "d1")

  # Initialize parameters
  add_args <- tryCatch(list(...), error = function(cond) list())
  init <- add_args$init
  init_params <- mfbvar:::parameter_initialization(Y = Y,
                                                   n_vars = n_vars,
                                                   n_lags = n_lags,
                                                   n_T_ = n_T_,
                                                   init = init,
                                                   n_determ = n_determ,
                                                   params)

  # Initialize storage
  mfbvar:::storage_initialization(init_params = init_params,
                                  params = params,
                                  envir = envir,
                                  n_vars = n_vars,
                                  n_lags = n_lags,
                                  n_reps = n_reps,
                                  n_thin = n_thin,
                                  n_T = n_T,
                                  n_T_ = n_T_,
                                  n_determ = n_determ,
                                  n_fcst = n_fcst)

  # Initialize ssng priors
  if (ssng) {
    init_ssng <- mfbvar:::ssng_initialization(prior_ng, s)
    mfbvar:::list_to_variables(init_ssng, envir, "c0", "c1", "s")
  } else {
    phi_mu <- matrix(0, 1, 1)
    lambda_mu <- matrix(0, 1, 1)
    omega <- matrix(diag(prior_psi_Omega), nrow = 1)
    c0 <- 0
    c1 <- 0
    s <- 0
  }

  roots <- vector("numeric", n_reps/n_thin)
  num_tries <- roots

  # For the posterior of Pi
  mcmc_ssng_iw(Y[-(1:n_lags),],Pi,Sigma,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,Omega_Pi,prior_Pi_mean,
               prior_S,D_mat,dt,d1,d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,n_reps,n_burnin,
               n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,ssng)
  if (verbose) {
    cat("\n")
  }

  return_obj <- mget(c(params, prior_params, retrieved_params))


  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}

mcmc_sampler.mfbvar_ss_iw <- function(x, ...) {

  out <- mfbvar_steadystate_iw(x, FALSE, ...)

  return(out)
}
mcmc_sampler.mfbvar_ssng_iw <- function(x, ...) {

  out <- mfbvar_steadystate_iw(x, TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_minn_iw <- function(x, ...){

  required_params <- c("Y", "n_lags", "n_burnin", "n_reps")
  prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                    "lambda1", "lambda3", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "n_thin", "Z_1")
  params <- c("Z", "Pi", "Sigma")

  check_required_params(x, "Y", "n_lags", "n_burnin", "n_reps")
  n_vars <- ncol(x$Y)

  prior_nu <- n_vars + 2
  priors <- prior_Pi_Sigma(lambda1 = x$lambda1, lambda2 = x$lambda3, prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
                           n_lags = x$n_lags, prior_nu = prior_nu)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  Y <- x$Y
  freq <- x$freq
  n_fcst <- x$n_fcst
  verbose <- x$verbose
  n_lags <- x$n_lags
  lambda4 <- x$lambda4

  # Add terms for constant
  prior_Pi_Omega <- diag(c(x$lambda1^2*lambda4^2, diag(prior_Pi_Omega)))
  prior_Pi_mean <- rbind(0, prior_Pi_mean)

  add_args <- list(...)
  n_reps <- x$n_reps
  n_burnin <- x$n_burnin
  n_thin <- ifelse(is.null(x$n_thin), 1, x$n_thin)
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_Z <- init$init_Z

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  freqs <- x$freqs
  Lambda_ <- x$Lambda_
  n_q <- sum(freq == freqs[1])
  if (n_q < n_vars) {
    T_b <- max(which(!apply(apply(Y[, freq == freqs[2], drop = FALSE], 2, is.na), 1, any)))
  } else {
    T_b <- nrow(Y)
  }
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
  }

  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags
  d <- matrix(1, nrow = nrow(Y), ncol = 1)
  post_nu <- n_T_ + prior_nu

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
    Z[,, 1] <- fill_na(Y)
  } else {
    if (all(dim(Z[,, 1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,, 1]), collapse = " x ")))
    }

  }

  ols_results <- ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = 1)

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
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  mcmc_minn_iw(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,
               Omega_Pi,prior_Pi_mean,prior_S,Z_1,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,
               n_thin,verbose,2)
  if (verbose) {
    cat("\n")
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = NULL, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, aggregation = x$aggregation, n_determ = 1,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = prior_nu, post_nu = prior_nu + n_T_, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = NULL, prior_psi_mean = NULL, n_reps = n_reps, n_burnin = n_burnin, n_thin = n_thin, Lambda_ = Lambda_, freq = freq,
                     init = list(init_Pi = Pi[,, n_reps/n_thin], init_Sigma = Sigma[,, n_reps/n_thin], init_Z = Z[,, n_reps/n_thin]))

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}


