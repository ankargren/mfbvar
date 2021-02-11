mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){

  check_required_params(x, "Y", "n_lags", "n_burnin", "n_reps")
  n_vars <- ncol(x$Y)

  prior_Pi_Omega <- create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  Y <- x$Y
  freq <- x$freq
  verbose <- x$verbose

  n_vars <- ncol(Y)
  n_lags <- x$n_lags
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst

  ## Priors

  priormu <- x$priormu
  priorphiidi <- x$priorphiidi
  priorphifac <- x$priorphifac
  priorsigmaidi <- x$priorsigmaidi
  priorsigmafac <- x$priorsigmafac
  priorfacload <- x$priorfacload
  restrict <- x$restrict

  if (length(priorsigmaidi) == 1) {
    priorsigmaidi <- rep(priorsigmaidi, n_vars)
  }
  if (length(priorsigmafac) == 1) {
    priorsigmafac <- rep(priorsigmafac, n_fac)
  }

  bmu <- priormu[1]
  Bmu <- priormu[2]^2

  Bsigma <- c(priorsigmaidi, priorsigmafac)

  B011inv <- 1/10^8
  B022inv <- 1/10^12

  armatau2 <- matrix(priorfacload^2, n_vars, n_fac) # priorfacload is scalar, or matrix

  armarestr <- matrix(FALSE, nrow = n_vars, ncol = n_fac)
  if (restrict == "upper") armarestr[upper.tri(armarestr)] <- TRUE
  armarestr <- matrix(as.integer(!armarestr), nrow = nrow(armarestr), ncol = ncol(armarestr)) # restrinv

  a0idi <- priorphiidi[1]
  b0idi <- priorphiidi[2]
  a0fac <- priorphifac[1]
  b0fac <- priorphifac[2]

  priorh0 <- rep(-1.0, n_vars + n_fac)

  ## Initials

  add_args <- list(...)
  n_reps <- x$n_reps
  n_burnin <- x$n_burnin
  n_thin <- ifelse(is.null(x$n_thin), 1, x$n_thin)

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

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- fill_na(Y)
  } else {
    init_Z <- init$init_Z
  }

  ### SV regressions
  if (is.null(init$init_mu)) {
    init_mu <- log(error_variance)
  } else {
    init_mu <- init$init_mu
  }
  if (is.null(init$init_sigma)) {
    init_sigma <- rep(0.75, n_vars + n_fac)
  } else {
    init_sigma <- init$init_sigma
  }
  if (is.null(init$init_phi)) {
    init_phi <- rep(0.2, n_vars + n_fac)
  } else {
    init_phi <- init$init_phi
  }

  ### Factors and loadings
  if (is.null(init$init_facload)) {
    init_facload <-  matrix(rnorm(n_vars*n_fac, sd = 0.5)^2, nrow=n_vars, ncol=n_fac)
  } else {
    init_facload <- init$init_facload
  }
  if (is.null(init$init_f)) {
    init_f <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_f <- init$init_f
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  t(cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = n_T_, ncol = n_vars+n_fac, byrow = TRUE)))
  } else {
    init_latent <- init$init_latent
  }
  if (is.null(init$init_latent0)) {
    init_latent0 <-  numeric(n_vars + n_fac)
  } else {
    init_latent0 <- init$init_latent0
  }

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

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags + 1, n_reps/n_thin))
  Z <- array(init_Z, dim = c(n_T, n_vars, n_reps/n_thin))
  Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }

  mu <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac),
                   dim = c(n_vars, n_fac, n_reps/n_thin))
  f <- array(matrix(init_f, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]

  aux <- matrix(0, 1, 1)
  global <- c(0)
  local <- matrix(0, 1, 1)
  a <- -1
  slice <- c(0)
  gig <- TRUE

  mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,h,
                         aux,global,local,slice,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                          a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                          armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,
                          n_vars,n_T_,n_fcst,n_thin,verbose,a,gig)
  if (verbose) {
    cat("\n")
  }
  return_obj <- list(Pi = Pi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     Lambda_ = Lambda_, aggregation = x$aggregation, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, Y = Y, Z_1 = Z_1, bmu = bmu,
                     Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
                     b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
                     B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
                     armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps, n_burnin = n_burnin,
                     n_q = n_q, T_b_ = T_b-n_lags, n_lags = n_lags,
                     n_vars = n_vars, n_T_ = n_T_, n_fcst = n_fcst,
                     n_thin = n_thin, verbose = verbose,
                     init = list(init_Pi = Pi[,, n_reps/n_thin],
                                 init_Z = Z[,, n_reps/n_thin],
                                 init_mu = mu[, n_reps/n_thin],
                                 init_phi = phi[, n_reps/n_thin],
                                 init_sigma = sigma[, n_reps/n_thin],
                                 init_facload = facload[,,n_reps/n_thin],
                                 init_f = f[,,n_reps/n_thin],
                                 init_h = h[,,n_reps/n_thin]))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)

}

mcmc_sampler.mfbvar_dl_fsv <- function(x, ...){

  check_required_params(x, "Y", "n_lags", "n_burnin", "n_reps")
  n_vars <- ncol(x$Y)

  prior_Pi_Omega <- create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  Y <- x$Y
  freq <- x$freq
  verbose <- x$verbose

  n_vars <- ncol(Y)
  n_lags <- x$n_lags
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst

  ## Priors

  priormu <- x$priormu
  priorphiidi <- x$priorphiidi
  priorphifac <- x$priorphifac
  priorsigmaidi <- x$priorsigmaidi
  priorsigmafac <- x$priorsigmafac
  priorfacload <- x$priorfacload
  restrict <- x$restrict

  if (length(priorsigmaidi) == 1) {
    priorsigmaidi <- rep(priorsigmaidi, n_vars)
  }
  if (length(priorsigmafac) == 1) {
    priorsigmafac <- rep(priorsigmafac, n_fac)
  }

  bmu <- priormu[1]
  Bmu <- priormu[2]^2

  Bsigma <- c(priorsigmaidi, priorsigmafac)

  B011inv <- 1/10^8
  B022inv <- 1/10^12

  armatau2 <- matrix(priorfacload^2, n_vars, n_fac) # priorfacload is scalar, or matrix

  armarestr <- matrix(FALSE, nrow = n_vars, ncol = n_fac)
  if (restrict == "upper") armarestr[upper.tri(armarestr)] <- TRUE
  armarestr <- matrix(as.integer(!armarestr), nrow = nrow(armarestr), ncol = ncol(armarestr)) # restrinv

  a0idi <- priorphiidi[1]
  b0idi <- priorphiidi[2]
  a0fac <- priorphifac[1]
  b0fac <- priorphifac[2]

  priorh0 <- rep(-1.0, n_vars + n_fac)

  ## DL
  if (!("a" %in% names(x))) {
    a <- 1
  } else {
    a <- x$a
  }

  gig <- ifelse(is.null(x$gig), TRUE, FALSE)

  RcppParallel::setThreadOptions(numThreads = x$n_cores)
  ## Initials

  add_args <- list(...)
  n_reps <- x$n_reps
  n_burnin <- x$n_burnin
  n_thin <- ifelse(is.null(x$n_thin), 1, x$n_thin)

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

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- fill_na(Y)
  } else {
    init_Z <- init$init_Z
  }

  ### SV regressions
  if (is.null(init$init_mu)) {
    init_mu <- log(error_variance)
  } else {
    init_mu <- init$init_mu
  }
  if (is.null(init$init_sigma)) {
    init_sigma <- rep(0.75, n_vars + n_fac)
  } else {
    init_sigma <- init$init_sigma
  }
  if (is.null(init$init_phi)) {
    init_phi <- rep(0.2, n_vars + n_fac)
  } else {
    init_phi <- init$init_phi
  }

  ### Factors and loadings
  if (is.null(init$init_facload)) {
    init_facload <-  matrix(rnorm(n_vars*n_fac, sd = 0.5)^2, nrow=n_vars, ncol=n_fac)
  } else {
    init_facload <- init$init_facload
  }
  if (is.null(init$init_f)) {
    init_f <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_f <- init$init_f
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  t(cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = n_T_, ncol = n_vars+n_fac, byrow = TRUE)))
  } else {
    init_latent <- init$init_latent
  }
  if (is.null(init$init_latent0)) {
    init_latent0 <-  numeric(n_vars + n_fac)
  } else {
    init_latent0 <- init$init_latent0
  }

  if (is.null(init$init_global)) {
    init_global <- 0.1
  } else {
    init_global <- init$init_global
  }

  if (is.null(init$init_aux)) {
    init_aux <- c(sqrt(prior_Pi_Omega[-1,])/init_global)
  } else {
    init_aux <- init$init_aux
  }

  if (is.null(init$init_local)) {
    init_local <- c(sqrt(prior_Pi_Omega[-1,])/init_global)
  } else {
    init_local <- init$init_local
  }

  if (is.null(init$init_slice)) {
    init_slice <- rep(1, n_vars^2*n_lags)
  } else {
    init_slice <- init$init_slice
  }

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

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags + 1, n_reps/n_thin))
  Z <- array(init_Z, dim = c(n_T, n_vars, n_reps/n_thin))
  Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }

  mu <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac),
                   dim = c(n_vars, n_fac, n_reps/n_thin))
  f <- array(matrix(init_f, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))
  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
             dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  aux <- matrix(init_aux, nrow = n_reps/n_thin, ncol = n_vars*n_vars*n_lags, byrow = TRUE)
  local <- matrix(init_local, nrow = n_reps/n_thin, ncol = n_vars*n_vars*n_lags, byrow = TRUE)
  global <- matrix(init_global, n_reps/n_thin, ncol = 1)
  slice <- matrix(init_slice, nrow = 1, ncol = n_vars*n_vars*n_lags)

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]
  mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,h,
                         aux,global,local,slice,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                         a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                         armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,
                         n_vars,n_T_,n_fcst,n_thin,verbose,a,gig)
  if (verbose) {
    cat("\n")
  }
  return_obj <- list(Pi = Pi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     aux = aux, local = local, global = global,
                     Lambda_ = Lambda_, aggregation = x$aggregation, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, Y = Y, Z_1 = Z_1, bmu = bmu,
                     Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
                     b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
                     B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
                     armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps, n_burnin = n_burnin,
                     n_q = n_q, T_b_ = T_b-n_lags, n_lags = n_lags,
                     n_vars = n_vars, n_T_ = n_T_, n_fcst = n_fcst,
                     n_thin = n_thin, verbose = verbose)

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)

}

mfbvar_steadystate_fsv <- function(x, ssng, ...) {

  if (ssng) {
    required_params <- c("Y", "d", "prior_psi_mean", "n_lags",
                          "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "prior_psi_mean",
                      "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                      "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict", "prior_ng", "s",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                      "n_T", "n_T_", "n_thin", "n_determ")
    params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")
  } else {
    required_params <- c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags",
                        "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "prior_psi_mean", "prior_psi_Omega",
                      "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                      "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                      "n_T", "n_T_", "n_thin", "n_determ")
    params <- c("Z", "psi", "Pi", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")
  }


  # Check inputs
  mfbvar:::check_required_params(x, required_params)

  # Assign variables
  mfbvar:::list_to_variables(x, parent.frame(), prior_params)

  # Retrieve some additional variables
  init_vars <- mfbvar:::variable_initialization(Y = Y, freq = freq, freqs = freqs,
                                       n_lags = n_lags, Lambda_ = Lambda_,
                                       n_thin = n_thin, d = d, d_fcst = d_fcst)
  mfbvar:::list_to_variables(init_vars, parent.frame(),
                    retrieved_params)

  # Prior
  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(lambda1, lambda2, lambda3, prior_Pi_AR1, Y, n_lags)[-1, ]

  # Initialize fsv priors
  init_fsv <- mfbvar:::fsv_initialization(priorsigmaidi = priorsigmaidi,
                                 priorsigmafac = priorsigmafac, priormu= priormu,
                                 restrict = restrict, priorphiidi = priorphiidi,
                                 priorphifac = priorphifac, n_vars = n_vars,
                                 n_fac = n_fac)
  mfbvar:::list_to_variables(init_fsv, parent.frame(), "priorsigmaidi", "priorsigmafac",
                    "bmu", "Bmu", "Bsigma", "B011inv", "B022inv", "armatau2",
                    "armarestr", "a0idi", "b0idi", "a0fac", "b0fac", "priorh0")

  # Initialize ssng priors
  if (ssng) {
    init_ssng <- mfbvar:::ssng_initialization(prior_ng, s)
    mfbvar:::list_to_variables(init_ssng, parent.frame(), "c0", "c1", "s")
  } else {
    phi_mu <- matrix(0, 1, 1)
    lambda_mu <- matrix(0, 1, 1)
    omega <- matrix(diag(prior_psi_Omega), nrow = 1)
    c0 <- 0
    c1 <- 0
    s <- 0
  }

  # Initialize parameters
  add_args <- list(...)
  init <- add_args$init
  init_params <- mfbvar:::parameter_initialization(Y = Y, n_vars = n_vars, n_lags = n_lags, n_T_ = n_T_,
                  init = init, n_fac = n_fac, n_determ = n_determ, params)

  # Initialize storage
  storage_initialization(init_params = init_params, params = params, envir = parent.frame(),
                         n_vars = n_vars, n_reps = n_reps, n_thin = n_thin,
                         n_T = n_T, n_T_ = n_T_, n_determ = n_determ,
                         n_fac = n_fac)


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
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]


  mcmc_ssng_fsv(Y[-(1:n_lags),],Pi,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,
                       mu,phi,sigma,f,facload,h,
                       Lambda_,prior_Pi_Omega,prior_Pi_AR1,D_mat,dt,d1,
                       d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,bmu,Bmu,
                       a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                       armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,
                       n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,TRUE)
  if (verbose) {
    cat("\n")
  }

  return_obj <- list(Pi = Pi, psi = psi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
              sigma = sigma, f = f, facload = facload, h = h,
              Lambda_ = Lambda_, aggregation = x$aggregation, prior_Pi_Omega = prior_Pi_Omega,
              prior_Pi_AR1 = prior_Pi_AR1, prior_psi_mean = prior_psi_mean,
              prior_psi_Omega = diag(omega[1, ]), d = d, Y = Y, Z_1 = Z_1, bmu = bmu,
              Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
              b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
              B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
              armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps,
              n_q = n_q, T_b_ = T_b-n_lags, n_lags = n_lags,
              n_vars = n_vars, n_T_ = n_T_, n_fcst = n_fcst, n_determ = n_determ,
              n_thin = n_thin, verbose = verbose,
              init = list(init_Pi = Pi[,, n_reps/n_thin],
                          init_psi = psi[n_reps/n_thin, ],
                          init_Z = Z[,, n_reps/n_thin],
                          init_mu = mu[, n_reps/n_thin],
                          init_phi = phi[, n_reps/n_thin],
                          init_sigma = sigma[, n_reps/n_thin],
                          init_facload = facload[,,n_reps/n_thin],
                          init_f = f[,,n_reps/n_thin],
                          init_h = h[,,n_reps/n_thin]))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)
}
mcmc_sampler.mfbvar_ss_fsv <- function(x, ...){

  out <- mfbvar_steadystate_fsv(x, FALSE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ssng_fsv <- function(x, ...){

  out <- mfbvar_steadystate_fsv(x, TRUE, ...)

  return(out)
}
