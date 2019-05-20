mcmc_sampler.mfbvar_minn_csv <- function(x, ...){
  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_nu <- n_vars + 2
  priors <- mfbvar:::prior_Pi_Sigma(lambda1 = x$lambda1, lambda2 = x$lambda3, prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
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

  phi_invvar <- 1/x$prior_phi[2]
  phi_meaninvvar <- x$prior_phi[1] * phi_invvar
  prior_sigma2 <- x$prior_sigma2[1]
  prior_df <- x$prior_sigma2[2]

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_Z <- init$init_Z
  init_phi <- init$init_phi
  init_sigma <- init$init_sigma
  init_f <- init$init_f

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
  phi <- rep(NA, n_reps/n_thin)
  sigma <- rep(NA, n_reps/n_thin)
  f <- matrix(NA, n_reps, n_T_)
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

  if (is.null(init_phi)) {
    phi[1] <- x$prior_phi[1]
  } else {
    phi[1] <- init_phi
  }

  if (is.null(init_sigma)) {
    sigma[1] <- sqrt(x$prior_sigma2[1])
  } else {
    sigma[1] <- init_sigma
  }

  if (is.null(init_f)) {
    f[1,] <- 0.0
  } else {
    f[1,] <- init_f
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]

  # For the posterior of Pi
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  set.seed(1)
  mfbvar:::mcmc_minn_csv(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,phi,sigma,f,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,
                        Omega_Pi,prior_Pi_mean,prior_S,Z_1,10,phi_invvar,phi_meaninvvar,prior_sigma2,prior_df,
                        n_reps,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,n_thin,verbose)

  return_obj <- list(Pi = Pi, Sigma = Sigma, Z = Z, phi = phi, sigma = sigma, f = f,
                     Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars,
                     n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_mean = prior_Pi_mean, prior_S = prior_S,
                     prior_nu = n_vars+2, post_nu = n_T + n_vars+2, d = d, Y = Y,
                     n_T = n_T, n_T_ = n_T_, n_reps = n_reps, Lambda_ = Lambda_,
                     init = list(init_Pi = Pi[,, n_reps/n_thin],
                                 init_Sigma = Sigma[,, n_reps/n_thin],
                                 init_Z = Z[,, n_reps/n_thin],
                                 init_phi = phi[n_reps/n_thin],
                                 init_sigma = sigma[n_reps/n_thin],
                                 init_f = f[n_reps/n_thin,]))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)

}

mcmc_sampler.mfbvar_ss_csv <- function(x, ...) {

  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$prior_psi_Omega) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }
  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }

  priors <- mfbvar:::prior_Pi_Sigma(lambda1 = x$lambda1, lambda2 = x$lambda3, prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
                                    n_lags = x$n_lags, prior_nu = n_vars + 2)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  Y <- x$Y
  d <- x$d
  d_fcst <- x$d_fcst
  freq <- x$freq
  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  n_fcst <- x$n_fcst
  check_roots <- x$check_roots
  verbose <- x$verbose

  phi_invvar <- 1/x$prior_phi[2]
  phi_meaninvvar <- x$prior_phi[1] * phi_invvar
  prior_sigma2 <- x$prior_sigma2[1]
  prior_df <- x$prior_sigma2[2]

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(is.null(add_args$n_thin),1,add_args$n_thin)
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_psi <- init$init_psi
  init_Z <- init$init_Z
  init_phi <- init$init_phi
  init_sigma <- init$init_sigma
  init_f <- init$init_f

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)
  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_Pi_mean)))/n_vars^2
  n_q <- sum(freq == "q")
  n_m <- sum(freq == "m")
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
  }
  y_in_p <- Y[-(1:n_lags), ]
  if (n_q < n_vars) {
    T_b <- min(apply(y_in_p[,1:n_m], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  } else {
    T_b <- nrow(y_in_p)
  }
  if (n_q > 0) {
    Lambda_ <- mfbvar:::build_Lambda(rep("q", n_q), 3)
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }

  n_pseudolags <- max(c(n_lags, 3))
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

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

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps/n_thin))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps/n_thin))
  psi   <- array(NA, dim = c(n_reps/n_thin, n_vars * n_determ))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps/n_thin))
  phi <- rep(NA, n_reps/n_thin)
  sigma <- rep(NA, n_reps/n_thin)
  f <- matrix(NA, n_reps/n_thin, n_T_)
  Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]
  roots <- vector("numeric", n_reps/n_thin)
  num_tries <- roots



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

  ols_results <- tryCatch(mfbvar:::ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ),
                          error = function(cond) NULL)
  if (is.null(ols_results)) {
    ols_results <- list()
    ols_results$Pi <- prior_Pi_mean
    ols_results$S <- prior_S
    ols_results$psi <- prior_psi_mean
  }

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- mfbvar:::build_companion(Pi = Pi[,, 1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- mfbvar:::max_eig_cpp(Pi_comp)
  }

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    if (all(dim(Sigma[,,1]) == dim(init_Sigma))) {
      Sigma[,, 1] <- init_Sigma
    } else {
      stop(paste0("The dimension of init_Sigma is ", paste(dim(init_Sigma), collapse = " x "), ", but should be ", paste(dim(Sigma[,,1]), collapse = " x ")))
    }
  }

  if (is.null(init_psi)) {
    if (roots[1] < 1) {
      psi[1, ] <- ols_results$psi
    } else {
      psi[1, ] <- prior_psi_mean
    }
  } else {
    if (length(psi[1, ]) == length(init_psi)) {
      psi[1,] <- init_psi
    } else {
      stop(paste0("The length of init_psi is ", paste(length(init_psi), collapse = " x "), ", but should be ", paste(length(psi[1,]), collapse = " x ")))
    }
  }

  if (is.null(init_phi)) {
    phi[1] <- x$prior_phi[1]
  } else {
    phi[1] <- init_phi
  }

  if (is.null(init_sigma)) {
    sigma[1] <- sqrt(x$prior_sigma2[1])
  } else {
    sigma[1] <- init_sigma
  }

  if (is.null(init_f)) {
    f[1,] <- 0.0
  } else {
    f[1,] <- init_f
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  # if requested
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]

  # For the posterior of Pi
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  # For the posterior of psi
  inv_prior_psi_Omega <- solve(prior_psi_Omega)
  inv_prior_psi_Omega_mean <- inv_prior_psi_Omega %*% prior_psi_mean
  Z_1 <- Z[1:n_pseudolags,, 1]

  mfbvar:::mcmc_ss_csv(Y[-(1:n_lags),],Pi,Sigma,psi,Z,Z_fcst,phi,sigma,f,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,Omega_Pi,prior_Pi_mean,
                       prior_S,D_mat,dt,d1,d_fcst_lags,inv_prior_psi_Omega,inv_prior_psi_Omega_mean,check_roots,Z_1,
                       10,phi_invvar,phi_meaninvvar,prior_sigma2,prior_df,n_reps,n_q,T_b,n_lags,n_vars,n_T_,n_fcst,n_determ,n_thin,verbose)

  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z, phi = phi, sigma = sigma, f = f, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, smoothed_Z = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = n_vars+2, post_nu = n_T + n_vars+2, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = prior_psi_Omega, prior_psi_mean = prior_psi_mean, n_reps = n_reps, Lambda_ = Lambda_,
                     init = list(init_Pi = Pi[,, n_reps/n_thin], init_Sigma = Sigma[,, n_reps/n_thin], init_psi = psi[n_reps/n_thin, ], init_Z = Z[,, n_reps/n_thin], init_phi = phi[n_reps/n_thin], init_sigma = sigma[n_reps/n_thin], init_f = f[n_reps/n_thin,]))

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)
}

mcmc_sampler.mfbvar_ssng_csv <- function(x, ...) {

  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }
  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }

  priors <- mfbvar:::prior_Pi_Sigma(lambda1 = x$lambda1, lambda2 = x$lambda3, prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
                                    n_lags = x$n_lags, prior_nu = n_vars + 2)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  Y <- x$Y
  d <- x$d
  d_fcst <- x$d_fcst
  freq <- x$freq
  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  n_fcst <- x$n_fcst
  check_roots <- x$check_roots
  verbose <- x$verbose

  phi_invvar <- 1/x$prior_phi[2]
  phi_meaninvvar <- x$prior_phi[1] * phi_invvar
  prior_sigma2 <- x$prior_sigma2[1]
  prior_df <- x$prior_sigma2[2]

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(is.null(add_args$n_thin),1,add_args$n_thin)
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_psi <- init$init_psi
  init_Z <- init$init_Z
  init_phi <- init$init_phi
  init_sigma <- init$init_sigma
  init_f <- init$init_f
  init_omega <- init$init_omega
  init_phi_mu <- init$init_phi_mu
  init_lambda_mu <- init$init_lambda_mu

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)
  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_Pi_mean)))/n_vars^2
  n_q <- sum(freq == "q")
  n_m <- sum(freq == "m")
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
  }
  y_in_p <- Y[-(1:n_lags), ]
  if (n_q < n_vars) {
    T_b <- min(apply(y_in_p[,1:n_m], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  } else {
    T_b <- nrow(y_in_p)
  }
  if (n_q > 0) {
    Lambda_ <- mfbvar:::build_Lambda(rep("q", n_q), 3)
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }
  n_pseudolags <- max(c(n_lags, 3))
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags



  c0 <- ifelse(is.null(x$c0), 0.01, x$c0)
  c1 <- ifelse(is.null(x$c1), 0.01, x$c1)
  s <- ifelse(is.null(x[["s"]]), -10, x$s)
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

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps/n_thin))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps/n_thin))
  psi   <- array(NA, dim = c(n_reps/n_thin, n_vars * n_determ))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps/n_thin))
  phi <- rep(NA, n_reps/n_thin)
  sigma <- rep(NA, n_reps/n_thin)
  f <- matrix(NA, n_reps/n_thin, n_T_)
  omega <- matrix(NA, nrow = n_reps/n_thin, ncol = n_vars * n_determ)
  phi_mu <- rep(NA, n_reps/n_thin)
  lambda_mu <- rep(NA, n_reps/n_thin)
  Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]
  roots <- vector("numeric", n_reps/n_thin)
  num_tries <- roots



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

  ols_results <- tryCatch(mfbvar:::ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ),
                          error = function(cond) NULL)
  if (is.null(ols_results)) {
    ols_results <- list()
    ols_results$Pi <- prior_Pi_mean
    ols_results$S <- prior_S
    ols_results$psi <- prior_psi_mean
  }

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- mfbvar:::build_companion(Pi = Pi[,, 1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- mfbvar:::max_eig_cpp(Pi_comp)
  }

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    if (all(dim(Sigma[,,1]) == dim(init_Sigma))) {
      Sigma[,, 1] <- init_Sigma
    } else {
      stop(paste0("The dimension of init_Sigma is ", paste(dim(init_Sigma), collapse = " x "), ", but should be ", paste(dim(Sigma[,,1]), collapse = " x ")))
    }
  }

  if (is.null(init_psi)) {
    if (roots[1] < 1) {
      psi[1, ] <- ols_results$psi
    } else {
      psi[1, ] <- prior_psi_mean
    }
  } else {
    if (length(psi[1, ]) == length(init_psi)) {
      psi[1,] <- init_psi
    } else {
      stop(paste0("The length of init_psi is ", paste(length(init_psi), collapse = " x "), ", but should be ", paste(length(psi[1,]), collapse = " x ")))
    }
  }

  if (is.null(init_phi)) {
    phi[1] <- x$prior_phi[1]
  } else {
    phi[1] <- init_phi
  }

  if (is.null(init_sigma)) {
    sigma[1] <- sqrt(x$prior_sigma2[1])
  } else {
    sigma[1] <- init_sigma
  }

  if (is.null(init_f)) {
    f[1,] <- 0.0
  } else {
    f[1,] <- init_f
  }

  if (is.null(init_omega)) {
    if (is.null(prior_psi_Omega)) {
      omega[1, ] <- diag(prior_psi_Omega)
    } else {
      omega[1, ] <- rep(0.1, n_determ*n_vars)
    }
  } else {
    omega[1, ] <- init_omega
  }

  if (is.null(init_phi_mu)) {
    phi_mu[1] <- 1
  } else {
    phi_mu[1] <- init_phi_mu
  }

  if (is.null(init_lambda_mu)) {
    lambda_mu[1] <- 1
  } else {
    lambda_mu[1] <- init_lambda_mu
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  # if requested
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]

  # For the posterior of Pi
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  Z_1 <- Z[1:n_pseudolags,, 1]

  mfbvar:::mcmc_ssng_csv(Y[-(1:n_lags),],Pi,Sigma,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,phi,sigma,f,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,Omega_Pi,prior_Pi_mean,
                       prior_S,D_mat,dt,d1,d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,
                       10,phi_invvar,phi_meaninvvar,prior_sigma2,prior_df,n_reps,n_q,T_b,n_lags,n_vars,n_T_,n_fcst,n_determ,n_thin,verbose)

  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z, phi_mu = phi_mu, lambda_mu = lambda_mu, omega = omega,
                     phi = phi, sigma = sigma, f = f, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, smoothed_Z = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = n_vars+2, post_nu = n_T + n_vars+2, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = prior_psi_Omega, prior_psi_mean = prior_psi_mean, n_reps = n_reps, Lambda_ = Lambda_,
                     init = list(init_Pi = Pi[,, n_reps/n_thin], init_Sigma = Sigma[,, n_reps/n_thin], init_psi = psi[n_reps/n_thin, ], init_Z = Z[,, n_reps/n_thin], init_phi = phi[n_reps/n_thin], init_sigma = sigma[n_reps/n_thin], init_f = f[n_reps/n_thin,], init_omega = omega[n_reps/n_thin, ], init_lambda_mu = lambda_mu[n_reps/n_thin], init_phi_mu = phi_mu[n_reps/n_thin]))

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)
}
