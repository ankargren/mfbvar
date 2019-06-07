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
  priorng <- x$priorng
  columnwise <- x$columnwise
  restrict <- x$restrict
  heteroskedastic <- x$heteroskedastic
  priorhomoskedastic <- x$priorhomoskedastic

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
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.005), n_fac, n_T_)
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

  latent <- array(init_latent, dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
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
                                                  quiet = TRUE, interweaving = 4, signswitch = TRUE,
                                                  startpara = startpara, startlatent = startlatent,
                                                  startlatent0 = startlatent0,
                                                  startfacload = startfacload, startfac = startfac, priormu = priormu,
                                                  priorphiidi = priorphiidi, priorphifac = priorphifac, priorsigmaidi = priorsigmaidi,
                                                  priorsigmafac = priorsigmafac, priorfacload = priorfacload, priorng = priorng,
                                                  columnwise = columnwise, restrict = restrict, heteroskedastic = heteroskedastic,
                                                  priorhomoskedastic = priorhomoskedastic), error = function(cond) cond)
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
                     facload = facload_storage, latent = latent,  mu = mu_storage, sigma = sigma_storage, phi = phi_storage,
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

#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_ss_fsv <- function(x, ...){

  if (is.null(x$n_fac)) {
    stop("The number of factors (n_fac) must be provided.")
  }

  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }

  Y <- x$Y
  d <- x$d
  d_fcst <- x$d_fcst
  n_vars <- ncol(Y)
  n_lags <- x$n_lags
  n_q <- sum(x$freq == "q")
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst
  n_determ <- dim(d)[2]
  mf <- TRUE
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
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

  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  check_roots <- x$check_roots
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  } else {
    num_tries <- NULL
  }


  init <- add_args$init
  error_variance <- mfbvar:::compute_error_variances(Y)

  priormu <- x$priormu
  priorphiidi <- x$priorphiidi
  priorphifac <- x$priorphifac
  priorsigmaidi <- x$priorsigmaidi
  priorsigmafac <- x$priorsigmafac
  priorfacload <- x$priorfacload
  priorng <- x$priorng
  columnwise <- x$columnwise
  restrict <- x$restrict
  heteroskedastic <- x$heteroskedastic
  priorhomoskedastic <- x$priorhomoskedastic

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*n_lags)
  } else {
    init_Pi <- init$init_Pi
  }


  ### Regression parameters
  if (is.null(init$init_psi)) {
    init_psi <- colMeans(y_in_p, na.rm = TRUE)
  } else {
    init_psi <- init$init_psi
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
    init_sigma <- rep(0.2, n_vars + n_fac)
  } else {
    init_sigma <- init$init_sigma
  }
  if (is.null(init$init_phi)) {
    init_phi <- rep(0.75, n_vars + n_fac)
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
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.005), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = n_T_, ncol = n_vars+n_fac, byrow = TRUE))
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

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, Y, n_lags)[-1, ]
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * n_T_) {
      par_fun <- mfbvar:::par_fun_top(rmvn_bcm)
    } else {
      par_fun <- mfbvar:::par_fun_top(rmvn_rue)
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

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags, n_reps/n_thin))
  Z <- array(init_Z, dim = c(n_T_, n_vars, n_reps/n_thin))
  psi   <- matrix(init_psi, n_reps, n_vars * n_determ, byrow = TRUE)
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                    dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]

  mu_storage <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/n_thin))
  fac_storage <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  latent <- array(init_latent, dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  Pi_i <- init_Pi
  Pi_i0 <- cbind(0, Pi_i)
  Z_i <- init_Z
  psi_i <- init_psi
  startpara <- list(mu = init_mu,
                    phi = init_phi,
                    sigma = init_sigma)
  startlatent <- latent[,,1]
  startlatent0 <- init_latent0
  startfacload <- matrix(init_facload, nrow = n_vars, ncol = n_fac)
  startfac <- matrix(init_fac, n_fac, n_T_)

  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]
  mu_mat <- dt %*% t(matrix(psi_i, nrow = n_vars))
  n_Lambda <- ncol(Lambda_)/nrow(Lambda_)
  mu_long <- matrix(0, n_Lambda+n_T_, n_vars)
  if (mf) {
    Lambda_single <- matrix(0, 1, n_Lambda)
    for (i in 1:n_Lambda) {
      Lambda_single[i] <- Lambda_[1, (i-1)*n_q+1]
    }
  }
  my <- matrix(0, nrow(y_in_p), ncol(y_in_p))

  if (verbose) {
    pb <- progress_bar$new(
      format = "[:bar] :percent eta: :eta",
      clear = FALSE, total = n_reps, width = 60)
  }

  error <- NULL

  inv_prior_psi_Omega <- solve(prior_psi_Omega)
  inv_prior_psi_Omega_mean <- inv_prior_psi_Omega %*% prior_psi_mean
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])


    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      mfbvar:::update_demean(my, mu_long, y_in_p, mu_mat, d1, matrix(psi_i, nrow = n_vars), Lambda_single, n_vars,
                    n_q, n_Lambda, n_T_)
    } else {
      mZ <- y_in_p - mu_mat
    }

    mZ1 <- Z_1 - d1 %*% t(matrix(psi_i, nrow = n_vars))
    Pi_i0[, -1] <- Pi_i

    if (mf){
      mZ <- tryCatch(mfbvar:::rsimsm_adaptive_univariate(my, Pi_i0, Sig, Lambda_, mZ1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i_demean, iter = i, block = "z")
        break
      }
    }
    Z_i_demean <- mZ
    Z_i <- Z_i_demean + mu_mat
    X <- mfbvar:::create_X_noint(rbind(Z_1, Z_i), n_lags)
    X_demean <- mfbvar:::create_X_noint(rbind(mZ1, Z_i_demean), n_lags)


    ## Produce forecasts


    ## Storage
    if (i %% n_thin == 0) {

      if (n_fcst > 0) {
        mu <- c(startpara$mu, numeric(n_fac))
        phi <- startpara$phi
        sigma <- startpara$sigma
        volatility_pred <- startlatent[n_T_, ]

        Z_pred <- matrix(0, n_fcst+n_lags, n_vars)
        Z_pred[1:n_lags, ] <- Z_i_demean[(n_T_-n_lags+1):n_T_,]
        for (j in 1:n_fcst) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- mfbvar:::create_X_t_noint(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
        Z_fcst[,,i/n_thin] <- Z_pred + d_fcst_lags %*% t(matrix(psi_i, nrow = n_vars))
      }

      Pi[,,i/n_thin] <- Pi_i
      Z[,,i/n_thin] <- Z_i
      psi[i/n_thin, ] <- psi_i

      mu_storage[,i/n_thin] <- startpara$mu
      sigma_storage[,i/n_thin] <- startpara$sigma
      phi_storage[,i/n_thin] <- startpara$phi

      fac_storage[,,i/n_thin] <- startfac
      facload_storage[,,i/n_thin] <- startfacload

      latent[,,i/n_thin] <- startlatent
    }

    ## Stochastic volatility block: sample latent factors, latent volatilities and factor loadings
    y_hat <- Z_i_demean - X_demean %*% t(Pi_i)
    fsample <- tryCatch(factorstochvol::fsvsample(y_hat, factors = n_fac, draws = 1, burnin = 0, priorh0idi = "stationary",
                                                  priorh0fac = "stationary", thin = 1, keeptime = "all",
                                                  runningstore = 0, runningstorethin = 10, runningstoremoments = 1,
                                                  quiet = TRUE, interweaving = 4, signswitch = TRUE,
                                                  startpara = startpara, startlatent = startlatent,
                                                  startlatent0 = startlatent0,
                                                  startfacload = startfacload, startfac = startfac, priormu = priormu,
                                                  priorphiidi = priorphiidi, priorphifac = priorphifac, priorsigmaidi = priorsigmaidi,
                                                  priorsigmafac = priorsigmafac, priorfacload = priorfacload, priorng = priorng,
                                                  columnwise = columnwise, restrict = restrict, heteroskedastic = heteroskedastic,
                                                  priorhomoskedastic = priorhomoskedastic), error = function(cond) cond)
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
    latent_nofac <- Z_i_demean - t(startfac) %*% t(startfacload)

    stationarity_check <- FALSE
    iter <- 0
    while(stationarity_check == FALSE) {
      iter <- iter + 1

      if (!parallelize) {
        if (prior_zero_mean) {
          for (j in 1:n_vars) {
            Pi_i[j,] <- tryCatch(mfbvar:::rmvn(X_demean/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5)), error = function(cond) cond)
          }
        } else {
          for (j in 1:n_vars) {
            Pi_i[j,] <- tryCatch(mfbvar:::rmvn_ccm(X_demean/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5), prior_Pi_AR1[j], j), error = function(cond) cond)
          }
        }
      } else {
        if (prior_zero_mean) {
          Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun, XX = X_demean, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac)), error = function(cond) cond)
        } else {
          Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun_AR1, XX = X_demean, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac, prior_Pi_AR1 = prior_Pi_AR1)), error = function(cond) cond)
        }
      }

      if (inherits(Pi_i, "error")) {
        warning("MCMC halted because of an error in the regression parameters step. See $error for more information.")
        error <- list(error = Pi_i, iter = i, block = "Pi_i")
        break
      }

      Pi_comp    <- mfbvar:::build_companion(Pi_i, n_vars = n_vars, n_lags = n_lags)
      if (check_roots == TRUE) {
        root <- mfbvar:::max_eig_cpp(Pi_comp)
      } else {
        root <- 0
      }
      if (root < 1) {
        stationarity_check <- TRUE
        if (check_roots == TRUE) {
          num_tries[i] <- iter
        }
      }
      if (iter == 1000) {
        warning("Attempted to draw stationary Pi 1,000 times.")
        error <- list(error = Pi_i, iter = i, block = "Pi_i")
        if (check_roots == TRUE) {
          num_tries[i] <- iter
        }
        break
      }

    }

    Pi_i0[, -1] <- Pi_i
    idivar <- exp(startlatent[, 1:n_vars])
    mfbvar:::posterior_psi_fsv(psi_i, mu_mat, Pi_i, D_mat, idivar, inv_prior_psi_Omega,
                  Z_i, X, startfacload, startfac, inv_prior_psi_Omega_mean, dt,
                  n_determ, n_vars, n_lags)

    if (verbose) {
      pb$tick()
    }
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Z = Z, psi = psi, Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst,
                     prior_Pi_Omega = prior_Pi_Omega, d = d, Y = Y, n_T = n_T, n_T_ = n_T_, n_reps = n_reps,
                     n_determ = n_determ, facload = facload_storage, latent = latent,  mu = mu_storage, sigma = sigma_storage, phi = phi_storage,
                     init = list(init_Pi = Pi_i, init_Z = Z_i, init_psi = psi_i, init_mu = startpara$mu,
                                 init_phi = startpara$phi, init_sigma = startpara$sigma,
                                 init_facload = startfacload,
                                 init_fac = startfac,
                                 init_latent = startlatent,
                                 init_latent0 = startlatent0),
                     num_tries = num_tries,
                     error = error)

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}

mcmc_sampler.mfbvar_ssng_fsv <- function(x, ...){

  if (is.null(x$n_fac)) {
    stop("The number of factors (n_fac) must be provided.")
  }

  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }

  Y <- x$Y
  d <- x$d
  d_fcst <- x$d_fcst
  n_vars <- ncol(Y)
  n_lags <- x$n_lags
  n_q <- sum(x$freq == "q")
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst
  n_determ <- dim(d)[2]
  mf <- TRUE
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    Y <- Y[complete_quarters, ]
    d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
    d <- d[complete_quarters, , drop = FALSE]
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

  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  check_roots <- x$check_roots
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  } else {
    num_tries <- NULL
  }
  c0 <- ifelse(is.null(x$c0), 0.01, x$c0)
  c1 <- ifelse(is.null(x$c1), 0.01, x$c1)
  s <- ifelse(is.null(x[["s"]]), -10, x$s)
  batch <- 0
  accept_vec <- numeric(n_reps)
  accept <- 0
  adaptive_mh <- FALSE
  if (s < 0) {
    M <- abs(s)
    s <- 1.0
    adaptive_mh <- TRUE
  }
  min_vec <- c(0.01, 0)

  init <- add_args$init
  error_variance <- mfbvar:::compute_error_variances(Y)

  priormu <- x$priormu
  priorphiidi <- x$priorphiidi
  priorphifac <- x$priorphifac
  priorsigmaidi <- x$priorsigmaidi
  priorsigmafac <- x$priorsigmafac
  priorfacload <- x$priorfacload
  priorng <- x$priorng
  columnwise <- x$columnwise
  restrict <- x$restrict
  heteroskedastic <- x$heteroskedastic
  priorhomoskedastic <- x$priorhomoskedastic

  ### Regression parameters
  if (is.null(init$init_Pi)) {
    init_Pi <- matrix(0, nrow = n_vars, ncol = n_vars*n_lags)
  } else {
    init_Pi <- init$init_Pi
  }


  ### Regression parameters
  if (is.null(init$init_psi)) {
    init_psi <- colMeans(y_in_p, na.rm = TRUE)
  } else {
    init_psi <- init$init_psi
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
    init_sigma <- rep(0.2, n_vars + n_fac)
  } else {
    init_sigma <- init$init_sigma
  }
  if (is.null(init$init_phi)) {
    init_phi <- rep(0.75, n_vars + n_fac)
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
    init_fac <-  matrix(rnorm(n_fac*n_T_, sd = 0.005), n_fac, n_T_)
  } else {
    init_fac <- init$init_fac
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = n_T_, ncol = n_vars+n_fac, byrow = TRUE))
  } else {
    init_latent <- init$init_latent
  }
  if (is.null(init$init_latent0)) {
    init_latent0 <-  numeric(n_vars + n_fac)
  } else {
    init_latent0 <- init$init_latent0
  }

  if (is.null(init$init_omega)) {
    if (is.null(prior_psi_Omega)) {
      init_omega <- diag(prior_psi_Omega)
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

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, Y, n_lags)[-1, ]
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * n_T_) {
      par_fun <- mfbvar:::par_fun_top(rmvn_bcm)
    } else {
      par_fun <- mfbvar:::par_fun_top(rmvn_rue)
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

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags, n_reps/n_thin))
  Z <- array(init_Z, dim = c(n_T_, n_vars, n_reps/n_thin))
  psi   <- matrix(init_psi, n_reps, n_vars * n_determ, byrow = TRUE)
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                    dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]

  omega <- matrix(init_omega, nrow = n_reps/n_thin, ncol = n_vars * n_determ)
  phi_mu <- rep(init_phi_mu, n_reps/n_thin)
  lambda_mu <- rep(init_lambda_mu, n_reps/n_thin)

  mu_storage <- matrix(init_mu, n_vars, n_reps/n_thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/n_thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/n_thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/n_thin))
  fac_storage <- array(matrix(init_fac, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin))

  latent <- array(init_latent, dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  Pi_i <- init_Pi
  Pi_i0 <- cbind(0, Pi_i)
  Z_i <- init_Z
  psi_i <- init_psi
  omega_i <- init_omega
  phi_mu_i <- init_phi_mu
  lambda_mu_i <- init_lambda_mu
  startpara <- list(mu = init_mu,
                    phi = init_phi,
                    sigma = init_sigma)
  startlatent <- latent[,,1]
  startlatent0 <- init_latent0
  startfacload <- matrix(init_facload, nrow = n_vars, ncol = n_fac)
  startfac <- matrix(init_fac, n_fac, n_T_)

  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]
  mu_mat <- dt %*% t(matrix(psi_i, nrow = n_vars))
  n_Lambda <- ncol(Lambda_)/nrow(Lambda_)
  mu_long <- matrix(0, n_Lambda+n_T_, n_vars)
  if (mf) {
    Lambda_single <- matrix(0, 1, n_Lambda)
    for (i in 1:n_Lambda) {
      Lambda_single[i] <- Lambda_[1, (i-1)*n_q+1]
    }
  }
  my <- matrix(0, nrow(y_in_p), ncol(y_in_p))

  if (verbose) {
    pb <- progress_bar$new(
      format = "[:bar] :percent eta: :eta",
      clear = FALSE, total = n_reps, width = 60)
  }

  error <- NULL
  inv_prior_psi_Omega <- solve(prior_psi_Omega)
  inv_prior_psi_Omega_mean <- inv_prior_psi_Omega %*% prior_psi_mean
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])

    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      mfbvar:::update_demean(my, mu_long, y_in_p, mu_mat, d1, matrix(psi_i, nrow = n_vars), Lambda_single, n_vars,
                             n_q, n_Lambda, n_T_)
    } else {
      mZ <- y_in_p - mu_mat
    }

    mZ1 <- Z_1 - d1 %*% t(matrix(psi_i, nrow = n_vars))
    Pi_i0[, -1] <- Pi_i

    if (mf){
      mZ <- tryCatch(mfbvar:::rsimsm_adaptive_univariate(my, Pi_i0, Sig, Lambda_, mZ1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i_demean, iter = i, block = "z")
        break
      }
    }
    Z_i_demean <- mZ
    Z_i <- Z_i_demean + mu_mat
    X <- mfbvar:::create_X_noint(rbind(Z_1, Z_i), n_lags)
    X_demean <- mfbvar:::create_X_noint(rbind(mZ1, Z_i_demean), n_lags)

    ## Produce forecasts


    ## Storage
    if (i %% n_thin == 0) {

      if (n_fcst > 0) {
        mu <- c(startpara$mu, numeric(n_fac))
        phi <- startpara$phi
        sigma <- startpara$sigma
        volatility_pred <- startlatent[n_T_, ]

        Z_pred <- matrix(0, n_fcst+n_lags, n_vars)
        Z_pred[1:n_lags, ] <- Z_i_demean[(n_T_-n_lags+1):n_T_,]
        for (j in 1:n_fcst) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- mfbvar:::create_X_t_noint(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
        Z_fcst[,,i/n_thin] <- Z_pred + d_fcst_lags %*% t(matrix(psi_i, nrow = n_vars))
      }

      Pi[,,i/n_thin] <- Pi_i
      Z[,,i/n_thin] <- Z_i
      psi[i/n_thin, ] <- psi_i

      mu_storage[,i/n_thin] <- startpara$mu
      sigma_storage[,i/n_thin] <- startpara$sigma
      phi_storage[,i/n_thin] <- startpara$phi

      fac_storage[,,i/n_thin] <- startfac
      facload_storage[,,i/n_thin] <- startfacload

      latent[,,i/n_thin] <- startlatent
    }

    ## Stochastic volatility block: sample latent factors, latent volatilities and factor loadings
    y_hat <- Z_i_demean - X_demean %*% t(Pi_i)
    fsample <- tryCatch(factorstochvol::fsvsample(y_hat, factors = n_fac, draws = 1, burnin = 0, priorh0idi = "stationary",
                                                  priorh0fac = "stationary", thin = 1, keeptime = "all",
                                                  runningstore = 0, runningstorethin = 10, runningstoremoments = 1,
                                                  quiet = TRUE, interweaving = 4, signswitch = TRUE,
                                                  startpara = startpara, startlatent = startlatent,
                                                  startlatent0 = startlatent0,
                                                  startfacload = startfacload, startfac = startfac, priormu = priormu,
                                                  priorphiidi = priorphiidi, priorphifac = priorphifac, priorsigmaidi = priorsigmaidi,
                                                  priorsigmafac = priorsigmafac, priorfacload = priorfacload, priorng = priorng,
                                                  columnwise = columnwise, restrict = restrict, heteroskedastic = heteroskedastic,
                                                  priorhomoskedastic = priorhomoskedastic), error = function(cond) cond)
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
    latent_nofac <- Z_i_demean - t(startfac) %*% t(startfacload)

    stationarity_check <- FALSE
    iter <- 0
    while(stationarity_check == FALSE) {
      iter <- iter + 1

      if (!parallelize) {
        if (prior_zero_mean) {
          for (j in 1:n_vars) {
            Pi_i[j,] <- tryCatch(mfbvar:::rmvn(X_demean/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5)), error = function(cond) cond)
          }
        } else {
          for (j in 1:n_vars) {
            Pi_i[j,] <- tryCatch(mfbvar:::rmvn_ccm(X_demean/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5), prior_Pi_AR1[j], j), error = function(cond) cond)
          }
        }
      } else {
        if (prior_zero_mean) {
          Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun, XX = X_demean, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac)), error = function(cond) cond)
        } else {
          Pi_i <- tryCatch(t(parallel::parSapply(cl, 1:n_vars, FUN = par_fun_AR1, XX = X_demean, startlatent = startlatent, D = prior_Pi_Omega, latent_nofac = latent_nofac, prior_Pi_AR1 = prior_Pi_AR1)), error = function(cond) cond)
        }
      }

      if (inherits(Pi_i, "error")) {
        warning("MCMC halted because of an error in the regression parameters step. See $error for more information.")
        error <- list(error = Pi_i, iter = i, block = "Pi_i")
        break
      }

      Pi_comp    <- mfbvar:::build_companion(Pi_i, n_vars = n_vars, n_lags = n_lags)
      if (check_roots == TRUE) {
        root <- mfbvar:::max_eig_cpp(Pi_comp)
      } else {
        root <- 0
      }
      if (root < 1) {
        stationarity_check <- TRUE
        if (check_roots == TRUE) {
          num_tries[i] <- iter
        }
      }
      if (iter == 1000) {
        warning("Attempted to draw stationary Pi 1,000 times.")
        error <- list(error = Pi_i, iter = i, block = "Pi_i")
        if (check_roots == TRUE) {
          num_tries[i] <- iter
        }
        break
      }

    }

    Pi_i0[, -1] <- Pi_i
    idivar <- exp(startlatent[, 1:n_vars])

    gig_lambda <- phi_mu_i-0.5
    gig_chi <- lambda_mu_i * phi_mu_i
    gig_psi <- (psi_i-prior_psi_mean)^2
    for (j in 1:(n_vars*n_determ)) {
      omega_i[j] = mfbvar:::do_rgig1(gig_lambda, gig_chi, gig_psi[j])
    }
    lambda_mu_i <- rgamma(1, n_vars*n_determ * phi_mu_i + c0, (0.5 * phi_mu_i * sum(omega_i) + c1))
    phi_mu_proposal <- phi_mu_i * exp(rnorm(1, sd = s))
    prob <- exp(mfbvar:::posterior_phi_mu(lambda_mu_i, phi_mu_proposal, omega_i, n_vars*n_determ)-mfbvar:::posterior_phi_mu(lambda_mu_i, phi_mu_i, omega_i, n_vars*n_determ)) * phi_mu_proposal/phi_mu_i
    u <- runif(1)
    if (u < prob) {
      phi_mu <- phi_mu_proposal
      accept <- 1
    } else {
      accept <- 0
    }
    if (adaptive_mh) {
      accept_vec[i] <- accept
      if (i %% 100 == 0) {
        batch <- batch + 1
        min_vec[2] <- batch^(-0.5)
        if (mean(accept_vec[(i-99):i]) > 0.44) {
          s_prop <- log(s) + min(min_vec)
          if (s_prop < M) {
            s <- exp(s_prop)
          }
        } else {
          s_prop <- log(s) - min(min_vec)
          if (s_prop > -M) {
            s <- exp(s_prop)
          }
        }
      }
    }


    mfbvar:::posterior_psi_fsv(psi_i, mu_mat, Pi_i, D_mat, idivar, inv_prior_psi_Omega,
                               Z_i, X, startfacload, startfac, inv_prior_psi_Omega_mean, dt,
                               n_determ, n_vars, n_lags)

    if (verbose) {
      pb$tick()
    }
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Z = Z, psi = psi, Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst,
                     prior_Pi_Omega = prior_Pi_Omega, d = d, Y = Y, n_T = n_T, n_T_ = n_T_, n_reps = n_reps,
                     n_determ = n_determ, facload = facload_storage, latent = latent,  mu = mu_storage, sigma = sigma_storage, phi = phi_storage,
                     init = list(init_Pi = Pi_i, init_Z = Z_i, init_psi = psi_i, init_mu = startpara$mu,
                                 init_phi = startpara$phi, init_sigma = startpara$sigma,
                                 init_facload = startfacload,
                                 init_fac = startfac,
                                 init_latent = startlatent,
                                 init_latent0 = startlatent0),
                     num_tries = num_tries,
                     error = error)

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}
