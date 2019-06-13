mcmc_sampler <- function(x, ...) {
  UseMethod("mcmc_sampler")
}


#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){
  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
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
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))

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
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
  }
  if (n_q > 0) {
    if (x$aggregation == "average") {
      Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)
    } else {
      Lambda_ <- mfbvar:::build_Lambda(rep("triangular", n_q), 5)}
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }


  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- mfbvar:::compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- mfbvar:::fill_na(Y)
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
  if (is.null(init$init_fac)) {
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
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
  f <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]

  aux <- matrix(0, 1, 1)
  global <- c(0)
  local <- matrix(0, 1, 1)
  a <- -1

  mfbvar:::mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,h,
                         aux,global,local,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                          a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                          armarestr,armatau2,n_fac,n_reps,n_q,T_b-n_lags,n_lags,
                          n_vars,n_T_,n_fcst,n_thin,verbose,a)

  return_obj <- list(Pi = Pi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     Lambda_ = Lambda_, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, Z_1 = Z_1, bmu = bmu,
                     Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
                     b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
                     B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
                     armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps,
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

#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_dl_fsv <- function(x, ...){
  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
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

  RcppParallel::setThreadOptions(numThreads = x$n_cores)

  ## Initials

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))

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
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
  }
  if (n_q > 0) {
    if (x$aggregation == "average") {
      Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)
    } else {
      Lambda_ <- mfbvar:::build_Lambda(rep("triangular", n_q), 5)}
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }


  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- mfbvar:::compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- mfbvar:::fill_na(Y)
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
  if (is.null(init$init_fac)) {
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
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
  f <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))
  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
             dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  aux <- matrix(init_aux, nrow = n_reps/n_thin, ncol = n_vars*n_vars*n_lags, byrow = TRUE)
  local <- matrix(init_local, nrow = n_reps/n_thin, ncol = n_vars*n_vars*n_lags, byrow = TRUE)
  global <- rep(init_global, n_reps/n_thin)

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]
  mfbvar:::mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,h,
                         aux,global,local,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                         a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                         armarestr,armatau2,n_fac,n_reps,n_q,T_b-n_lags,n_lags,
                         n_vars,n_T_,n_fcst,n_thin,verbose,a)

  return_obj <- list(Pi = Pi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     Lambda_ = Lambda_, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, Z_1 = Z_1, bmu = bmu,
                     Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
                     b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
                     B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
                     armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps,
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
                                 init_h = h[,,n_reps/n_thin],
                                 init_aux = aux,
                                 init_local = local,
                                 init_global = global))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)

}

mcmc_sampler.mfbvar_ss_fsv <- function(x, ...){
  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  Y <- x$Y
  freq <- x$freq
  verbose <- x$verbose

  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  d <- x$d
  d_fcst <- x$d_fcst
  check_roots <- x$check_roots
  n_determ <- dim(d)[2]

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
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))

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
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
  }
  if (n_q > 0) {
    if (x$aggregation == "average") {
      Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)
    } else {
      Lambda_ <- mfbvar:::build_Lambda(rep("triangular", n_q), 5)}
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }


  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- mfbvar:::compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Steady-states
  if (is.null(init$init_Z)) {
    init_Z <- mfbvar:::fill_na(Y)
  } else {
    init_Z <- init$init_Z
  }

  ### Latent high-frequency
  if (is.null(init$init_psi)) {
    init_psi <- prior_psi_mean
  } else {
    init_psi <- init$init_psi
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
  if (is.null(init$init_fac)) {
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
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
  psi <- array(init_psi, dim = c(n_reps/n_thin, n_vars * n_determ))
  Z <- array(init_Z, dim = c(n_T, n_vars, n_reps/n_thin))
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

  mu <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac),
                   dim = c(n_vars, n_fac, n_reps/n_thin))
  f <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
             dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]

  phi_mu <- matrix(0, 1, 1)
  lambda_mu <- matrix(0, 1, 1)
  omega <- matrix(diag(prior_psi_Omega), nrow = 1)
  c0 <- 0
  c1 <- 0
  s <- 0

  mfbvar:::mcmc_ss_fsv(Y[-(1:n_lags),],Pi,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,
                         mu,phi,sigma,f,facload,h,
                         Lambda_,prior_Pi_Omega,prior_Pi_AR1,D_mat,dt,d1,
                         d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,bmu,Bmu,
                         a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                         armarestr,armatau2,n_fac,n_reps,n_q,T_b-n_lags,n_lags,
                         n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,FALSE)

  return_obj <- list(Pi = Pi, psi = psi, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     Lambda_ = Lambda_, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, prior_psi_mean = prior_psi_mean,
                     prior_psi_Omega = diag(omega[1, ]), Z_1 = Z_1, bmu = bmu,
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

mcmc_sampler.mfbvar_ssng_fsv <- function(x, ...){
  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements: ", paste(names(test_sub)[which(test_sub)], collapse = " "))
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, x$n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  Y <- x$Y
  freq <- x$freq
  verbose <- x$verbose

  prior_psi_mean <- x$prior_psi_mean
  d <- x$d
  d_fcst <- x$d_fcst
  check_roots <- x$check_roots
  n_determ <- dim(d)[2]

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

  c0 <- ifelse(is.null(x$c0), 0.01, x$c0)
  c1 <- ifelse(is.null(x$c1), 0.01, x$c1)
  s <- ifelse(is.null(x[["s"]]), 1, x$s)

  ## Initials

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))

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
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
  }
  if (n_q > 0) {
    if (x$aggregation == "average") {
      Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)
    } else {
      Lambda_ <- mfbvar:::build_Lambda(rep("triangular", n_q), 5)}
  } else {
    Lambda_ <- matrix(0, 1, 3)
  }


  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  ## Initials
  init <- add_args$init
  y_in_p <- Y[-(1:n_lags), ]
  error_variance <- mfbvar:::compute_error_variances(Y)

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- mfbvar:::fill_na(Y)
  } else {
    init_Z <- init$init_Z
  }

  ### Steady-states
  if (is.null(init$init_psi)) {
    init_psi <- prior_psi_mean
  } else {
    init_psi <- init$init_psi
  }

  if (is.null(init$init_omega)) {
    if (!is.null(x$prior_psi_Omega)) {
      init_omega <- diag(x$prior_psi_Omega)
    } else {
      init_omega <- rep(0.1, n_determ*n_vars)
    }
  } else {
    init_omega <- init$init_omega
  }

  if (is.null(init$init_phi_mu)) {
    init_phi_mu <- 1
  } else {
    init_phi_mu <- init$init_phi_mu
  }

  if (is.null(init$init_lambda_mu)) {
    init_lambda_mu <- 1
  } else {
    init_lambda_mu <- init$init_lambda_mu
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
  if (is.null(init$init_fac)) {
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.5), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
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
  psi <- array(init_psi, dim = c(n_reps/n_thin, n_vars * n_determ))
  Z <- array(init_Z, dim = c(n_T, n_vars, n_reps/n_thin))
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

  mu <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac),
                   dim = c(n_vars, n_fac, n_reps/n_thin))
  f <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  h <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
             dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  omega <- matrix(init_omega, nrow = n_reps/n_thin, ncol = n_vars * n_determ, byrow = TRUE)
  phi_mu <- rep(init_phi_mu, n_reps/n_thin)
  lambda_mu <- rep(init_lambda_mu, n_reps/n_thin)

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]



  mfbvar:::mcmc_ss_fsv(Y[-(1:n_lags),],Pi,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,
                       mu,phi,sigma,f,facload,h,
                       Lambda_,prior_Pi_Omega,prior_Pi_AR1,D_mat,dt,d1,
                       d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,bmu,Bmu,
                       a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                       armarestr,armatau2,n_fac,n_reps,n_q,T_b-n_lags,n_lags,
                       n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,TRUE)

  return_obj <- list(Pi = Pi, psi = psi, omega = omega, lambda_mu = lambda_mu,
                     phi_mu = phi_mu, Z = Z, Z_fcst = NULL, mu = mu, phi = phi,
                     sigma = sigma, f = f, facload = facload, h = h,
                     Lambda_ = Lambda_, prior_Pi_Omega = prior_Pi_Omega,
                     prior_Pi_AR1 = prior_Pi_AR1, prior_psi_mean = prior_psi_mean,
                     prior_psi_Omega = diag(omega[1, ]), Z_1 = Z_1, bmu = bmu,
                     Bmu = Bmu, a0idi = a0idi, b0idi = b0idi, a0fac = a0fac,
                     b0fac = b0fac, Bsigma = Bsigma, B011inv = B011inv,
                     B022inv = B022inv, priorh0 = priorh0, armarestr = armarestr,
                     armatau2 = armatau2, n_fac = n_fac, n_reps = n_reps,
                     n_q = n_q, T_b_ = T_b-n_lags, n_lags = n_lags,
                     n_vars = n_vars, n_T_ = n_T_, n_fcst = n_fcst, n_determ = n_determ,
                     n_thin = n_thin, verbose = verbose,
                     init = list(init_Pi = Pi[,, n_reps/n_thin],
                                 init_psi = psi[n_reps/n_thin, ],
                                 init_omega = omega[n_reps/n_thin, ],
                                 init_lambda_mu = lambda_mu[n_reps/n_thin],
                                 init_phi_mu = phi_mu[n_reps/n_thin],
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
