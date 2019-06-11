#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){

  if (is.null(x$n_fac)) {
    stop("The number of factors (n_fac) must be provided.")
  }
  Y <- x$Y
  n_vars <- ncol(Y)
  n_lags <- x$n_lags
  n_q <- sum(x$freq == "q")
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst
  mf <- TRUE
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    mf <- FALSE
  }
  y_in_p <- Y[-(1:n_lags), ]
  if (n_q < n_vars) {
    T_b <- min(apply(y_in_p[,1:n_m], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  } else {
    T_b <- nrow(y_in_p)
  }
  n_T_ <- nrow(Y) - n_lags
  n_T <- nrow(Y)

  add_args <- list(...)
  n_reps <- add_args$n_reps
  n_thin <- ifelse(!is.null(add_args$n_thin), add_args$n_thin, ifelse(!is.null(x$n_thin), x$n_thin, 1))

  init <- add_args$init
  error_variance <- mfbvar:::compute_error_variances(Y)

  priormu <- x$priormu
  priorphiidi <- x$priorphiidi
  priorphifac <- x$priorphifac
  priorsigmaidi <- x$priorsigmaidi
  priorsigmafac <- x$priorsigmafac
  priorfacload <- x$priorfacload
  restrict <- x$restrict

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags + 1))
  } else {
    init_Pi <- init$init_Pi
  }

  ### Latent high-frequency
  if (is.null(init$init_Z)) {
    init_Z <- y_in_p
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
    init_facload <-  matrix(rnorm(n_vars*n_fac, sd = .5)^2, nrow=n_vars, ncol=n_fac)
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

  cl <- x$cl
  Z_1 <- mfbvar:::fill_na(Y)[(1:n_lags), ]
  verbose <- x$verbose

  ## Set up cluster (if used)
  if (!is.null(cl)) {
    parallelize <- TRUE
    parallel::clusterCall(cl, fun = function() library(mfbvar))
    parallel::clusterExport(cl, varlist = c("par_fun"))
  } else {
    parallelize <- FALSE
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, Y, n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * n_T_) {
      par_fun <- mfbvar:::par_fun_top(mfbvar:::rmvn_bcm)
    } else {
      par_fun <- mfbvar:::par_fun_top(mfbvar:::rmvn_rue)
    }
  }

  ## Obtain the aggregation matrix for the quarterly only
  if (mf) {
    if (x$aggregation == "average") {
      Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)
    } else {
      Lambda_ <- mfbvar:::build_Lambda(rep("triangular", n_q), 5)
    }
  }

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags + 1, n_reps/n_thin))
  Z <- array(init_Z, dim = c(n_T_, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                    dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }

  mu_storage <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/n_thin))
  fac_storage <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  latent <- array(t(init_latent), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  Pi_i <- Pi[,,1]
  Z_i <- init_Z
  startpara <- list(mu = init_mu,
                    phi = init_phi,
                    sigma = init_sigma)
  startlatent <- latent[,,1]
  startlatent0 <- init_latent0
  startfacload <- matrix(init_facload, nrow = n_vars, ncol = n_fac)
  startfac <- matrix(init_fac, n_fac, n_T_)

  if (verbose) {
    pb <- progress_bar$new(
      format = "[:bar] :percent eta: :eta",
      clear = FALSE, total = n_reps, width = 60)
  }

  if (!mf) {
    X <- mfbvar:::create_X(rbind(Z_1, Z_i), n_lags)
  }

  error <- NULL
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])

    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      Z_i <- tryCatch(mfbvar:::rsimsm_adaptive_univariate(y_in_p, Pi_i, Sig, Lambda_, Z_1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i, iter = i, block = "z")
        break
      }
      Z_i <- rbind(Z_1, Z_i)
      X <- mfbvar:::create_X(Z_i, n_lags)
      Z_i <- Z_i[-(1:n_lags), ]
    }

    ## Produce forecasts


    ## Storage
    if (i %% n_thin == 0) {

      if (n_fcst > 0) {
        mu <- c(startpara$mu, numeric(n_fac))
        phi <- startpara$phi
        sigma <- startpara$sigma
        volatility_pred <- startlatent[n_T_, ]

        Z_pred <- matrix(0, n_fcst+n_lags, n_vars)
        Z_pred[1:n_lags, ] <- Z_i[(n_T_-n_lags+1):n_T_,]
        for (j in 1:n_fcst) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- mfbvar:::create_X_t(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
        Z_fcst[,,i/n_thin] <- Z_pred
      }

      Pi[,,i/n_thin] <- Pi_i
      Z[,,i/n_thin] <- Z_i

      mu_storage[,i/n_thin] <- startpara$mu
      sigma_storage[,i/n_thin] <- startpara$sigma
      phi_storage[,i/n_thin] <- startpara$phi

      fac_storage[,,i/n_thin] <- startfac
      facload_storage[,,i/n_thin] <- startfacload

      latent[,,i/n_thin] <- startlatent
    }

    ## Stochastic volatility block: sample latent factors, latent volatilities and factor loadings
    y_hat <- Z_i - X %*% t(Pi_i)
    fsample <- tryCatch(factorstochvol::fsvsample(y_hat, factors = n_fac, draws = 1, burnin = 0, priorh0idi = "stationary",
                                                  priorh0fac = "stationary", thin = 1, keeptime = "all",
                                                  runningstore = 0, runningstorethin = 10, runningstoremoments = 1,
                                                  quiet = TRUE, interweaving = 4, signswitch = FALSE,
                                                  startpara = startpara, startlatent = startlatent,
                                                  startlatent0 = startlatent0,
                                                  startfacload = startfacload, startfac = startfac, priormu = priormu,
                                                  priorphiidi = priorphiidi, priorphifac = priorphifac, priorsigmaidi = priorsigmaidi,
                                                  priorsigmafac = priorsigmafac, priorfacload = priorfacload, restrict = restrict),
                                                  error = function(cond) cond)
    if (inherits(fsample, "error")) {
      warning("MCMC halted because of an error in the factor stochastic volatility step. See $error for more information.")
      error <- list(error = fsample, iter = i, block = "fsample")
      break
    }
    startpara$mu <- fsample$para[1,1:n_vars,1]
    startpara$phi <- fsample$para[2,,1]
    startpara$sigma <- fsample$para[3,,1]
    startlatent0 <- c(fsample$h0)
    startlatent <- fsample$h[,,1]
    startfacload <- matrix(fsample$facload[,,1], nrow = n_vars, ncol = n_fac)
    startfac <- matrix(fsample$f[,,1], nrow = n_fac)

    ## Regression parameters block: sample Pi (possibly in parallel)
    latent_nofac <- Z_i - t(startfac) %*% t(startfacload)

    if (!parallelize) {
      if (prior_zero_mean) {
        for (j in 1:n_vars) {
          Pi_i[j,] <- tryCatch(mfbvar:::rmvn(X/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5)), error = function(cond) cond)
        }
      } else {
        for (j in 1:n_vars) {
          Pi_i[j,] <- tryCatch(mfbvar:::rmvn_ccm(X/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5), prior_Pi_AR1[j], j), error = function(cond) cond)
        }
      }
    } else {
      if (prior_zero_mean) {
        Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun, XX = X, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac)), error = function(cond) cond)
      } else {
        Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun_AR1, XX = X, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac, prior_Pi_AR1 = prior_Pi_AR1)), error = function(cond) cond)
      }
    }
    if (inherits(Pi_i, "error")) {
      warning("MCMC halted because of an error in the regression parameters step. See $error for more information.")
      error <- list(error = Pi_i, iter = i, block = "Pi_i")
      break
    }

    if (verbose) {
      pb$tick()
    }
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Z = Z, Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst,
                     prior_Pi_Omega = prior_Pi_Omega, Y = Y, n_T = n_T, n_T_ = n_T_, n_reps = n_reps,
                     facload = facload_storage, f = fac_storage,
                     h = latent,  mu = mu_storage, sigma = sigma_storage, phi = phi_storage,
                     init = list(init_Pi = Pi_i, init_Z = Z_i, init_mu = startpara$mu,
                                 init_phi = startpara$phi, init_sigma = startpara$sigma,
                                 init_facload = startfacload,
                                 init_fac = startfac,
                                 init_latent = startlatent,
                                 init_latent0 = startlatent0),
                     error = error)

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}

mcmc_sampler.mfbvar_minn_fsv2 <- function(x, ...){
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

  mfbvar:::mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,h,
                          Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                          a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                          armarestr,armatau2,n_fac,n_reps,n_q,T_b-n_lags,n_lags,
                          n_vars,n_T_,n_fcst,n_thin,verbose)

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
