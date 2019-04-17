#' @rdname mcmc_sampler
mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){

  if (is.null(x$n_fac)) {
    stop("The number of factors (n_fac) must be provided.")
  }

  n_vars <- ncol(x$Y)
  n_lags <- x$n_lags
  n_q <- sum(x$freq == "q")
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst

  y_in_p <- x$Y[-(1:n_lags), ]

  T_b <- min(apply(y_in_p[,1:n_m], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  TT <- nrow(x$Y) - n_lags
  n_T <- nrow(x$Y)



  add_args <- list(...)

  n_reps <- add_args$n_reps
  if (!is.null(x$thin)) {
    thin <- x$thin
  } else {
    thin <- 1
  }

  init <- add_args$init
  error_variance <- mfbvar:::compute_error_variances(x$Y)

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
    init_fac <-  matrix(rnorm(n_fac*TT, sd = 0.005), n_fac, TT)
  } else {
    init_fac <- init$init_fac
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  t(cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = TT, ncol = n_vars+n_fac, byrow = TRUE)))
  } else {
    init_latent <- init$init_latent
  }
  if (is.null(init$init_latent0)) {
    init_latent0 <-  numeric(n_vars + n_fac)
  } else {
    init_latent0 <- init$init_latent0
  }

  cl <- x$cl
  Z1 <- mfbvar:::fill_na(x$Y)[(1:n_lags), ]
  verbose <- x$verbose
  mf <- TRUE

  ## Set up cluster (if used)
  if (!is.null(cl)) {
    parallelize <- TRUE
    parallel::clusterCall(cl, fun = function() library(mfbvar))
    parallel::clusterExport(cl, varlist = c("par_fun"))
  } else {
    parallelize <- FALSE
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, n_lags)
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * TT) {
      par_fun <- mfbvar:::par_fun_top(mfbvar:::rmvn_bcm)
    } else {
      par_fun <- mfbvar:::par_fun_top(mfbvar:::rmvn_rue)
    }
  }

  ## Obtain the aggregation matrix for the quarterly only
  if (mf) {
    Lambda_companion <- mfbvar:::build_Lambda(c(rep("m", n_vars-n_q), rep("q", n_q)), n_lags)
    Lambda_comp <- matrix(Lambda_companion[(n_m+1):n_vars, c(t(sapply((n_m+1):n_vars, function(x) seq(from = x, to = n_vars*n_lags, by = n_vars))))],
                          nrow = n_q)
  }

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags + 1, n_reps/thin))
  Z <- array(init_Z, dim = c(TT, n_vars, n_reps/thin))
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                    dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }

  mu_storage <- matrix(init_mu, n_vars, n_reps/thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/thin))
  fac_storage <- array(matrix(init_fac, n_fac, TT), dim = c(n_fac, TT, n_reps/thin))

  latent <- array(init_latent, dim = c(TT, n_vars+n_fac, n_reps/thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  Pi_i <- Pi[,,1]
  Z_i <- init_Z
  startpara <- list(mu = init_mu,
                    phi = init_phi,
                    sigma = init_sigma)
  startlatent <- latent[,,1]
  startlatent0 <- init_latent0
  startfacload <- matrix(init_facload, nrow = n_vars, ncol = n_fac)
  startfac <- matrix(init_fac, n_fac, TT)

  if (verbose) {
    pb <- progress_bar$new(
      format = "[:bar] :percent eta: :eta",
      clear = FALSE, total = n_reps, width = 60)
  }

  error <- NULL
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])

    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      Z_i <- tryCatch(mfbvar:::rsimsm_adaptive_univariate(y_in_p, Pi_i, Sig, Lambda_comp, Z1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i, iter = i, block = "z")
        break
      }
      Z_i <- rbind(Z1, Z_i)
      X <- mfbvar:::create_X(Z_i, n_lags)
      Z_i <- Z_i[-(1:n_lags), ]
    }

    ## Produce forecasts


    ## Storage
    if (i %% thin == 0) {

      if (n_fcst > 0) {
        mu <- c(startpara$mu, numeric(n_fac))
        phi <- startpara$phi
        sigma <- startpara$sigma
        volatility_pred <- startlatent[TT, ]

        Z_pred <- matrix(0, n_fcst+n_lags, n_vars)
        Z_pred[1:n_lags, ] <- Z_i[(TT-n_lags+1):TT,]
        for (j in 1:n_fcst) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- mfbvar:::create_X_t(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
        Z_fcst[,,i/thin] <- Z_pred
      }

      Pi[,,i/thin] <- Pi_i
      Z[,,i/thin] <- Z_i

      mu_storage[,i/thin] <- startpara$mu
      sigma_storage[,i/thin] <- startpara$sigma
      phi_storage[,i/thin] <- startpara$phi

      fac_storage[,,i/thin] <- startfac
      facload_storage[,,i/thin] <- startfacload

      latent[,,i/thin] <- startlatent
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
                     prior_Pi_Omega = prior_Pi_Omega, Y = x$Y, n_T = n_T, n_T_ = TT, n_reps = n_reps-1,
                     facload = facload_storage, latent = latent,  mu = mu, sigma = sigma, phi = phi,
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

  add_args <- list(...)
  n_reps <- add_args$n_reps

  if (is.null(x$n_fac)) {
    stop("The number of factors (n_fac) must be provided.")
  }

  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }


  d <- x$d
  d_fcst <- x$d_fcst
  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  check_roots <- x$check_roots
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  } else {
    num_tries <- NULL
  }

  n_vars <- ncol(x$Y)
  n_lags <- x$n_lags
  n_q <- sum(x$freq == "q")
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst
  n_determ <- dim(d)[2]

  y_in_p <- x$Y[-(1:n_lags), ]

  T_b <- min(apply(y_in_p[,1:n_m,drop=FALSE], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  TT <- nrow(x$Y) - n_lags
  n_T <- nrow(x$Y)






  if (!is.null(x$thin)) {
    thin <- x$thin
  } else {
    thin <- 1
  }

  init <- add_args$init
  error_variance <- mfbvar:::compute_error_variances(x$Y)

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
    init_fac <-  matrix(rnorm(n_fac*TT, sd = 0.005), n_fac, TT)
  } else {
    init_fac <- init$init_fac
  }

  ### Latent volatilities
  if (is.null(init$init_latent)) {
    init_latent <-  t(cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = TT, ncol = n_vars+n_fac, byrow = TRUE)))
  } else {
    init_latent <- init$init_latent
  }
  if (is.null(init$init_latent0)) {
    init_latent0 <-  numeric(n_vars + n_fac)
  } else {
    init_latent0 <- init$init_latent0
  }

  cl <- x$cl
  Z1 <- mfbvar:::fill_na(x$Y)[(1:n_lags), ]
  verbose <- x$verbose
  mf <- TRUE

  ## Set up cluster (if used)
  if (!is.null(cl)) {
    parallelize <- TRUE
    parallel::clusterCall(cl, fun = function() library(mfbvar))
    parallel::clusterExport(cl, varlist = c("par_fun"))
  } else {
    parallelize <- FALSE
  }

  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(x$lambda1, x$lambda2, x$lambda3, x$prior_Pi_AR1, x$Y, n_lags)[-1, ]
  prior_Pi_AR1 <- x$prior_Pi_AR1
  prior_zero_mean <- all(x$prior_Pi_AR1 == 0)

  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * TT) {
      par_fun <- mfbvar:::par_fun_top(rmvn_bcm)
    } else {
      par_fun <- mfbvar:::par_fun_top(rmvn_rue)
    }
  }

  ## Obtain the aggregation matrix for the quarterly only
  if (mf) {
    Lambda_companion <- mfbvar:::build_Lambda(c(rep("m", n_vars-n_q), rep("q", n_q)), n_lags)
    Lambda_comp <- matrix(Lambda_companion[(n_m+1):n_vars, c(t(sapply((n_m+1):n_vars, function(x) seq(from = x, to = n_vars*n_lags, by = n_vars))))],
                          nrow = n_q)
  }

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags, n_reps/thin))
  Z <- array(init_Z, dim = c(TT, n_vars, n_reps/thin))
  psi   <- matrix(init_psi, n_reps, n_vars * n_determ, byrow = TRUE)
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                    dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }

  mu_storage <- matrix(init_mu, n_vars, n_reps/thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/thin))
  fac_storage <- array(matrix(init_fac, n_fac, TT), dim = c(n_fac, TT, n_reps/thin))

  latent <- array(init_latent, dim = c(TT, n_vars+n_fac, n_reps/thin),
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
  startfac <- matrix(init_fac, n_fac, TT)

  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]
  mu_mat <- dt %*% t(matrix(psi_i, nrow = n_vars))

  if (verbose) {
    pb <- progress_bar$new(
      format = "[:bar] :percent eta: :eta",
      clear = FALSE, total = n_reps, width = 60)
  }

  error <- NULL
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])

    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      mZ <- y_in_p - mu_mat
      mZ <- as.matrix(mZ)
      mZ1 <- Z1 - d1 %*% t(matrix(psi_i, nrow = n_vars))
      Z_i_demean <- tryCatch(mfbvar:::rsimsm_adaptive_univariate(mZ, Pi_i0, Sig, Lambda_comp, mZ1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i_demean, iter = i, block = "z")
        break
      }
      Z_i <- Z_i_demean + mu_mat
      X <- mfbvar:::create_X_noint(rbind(Z1, Z_i), n_lags)
      #Z_i_demean <- Z_i
      X_demean <- mfbvar:::create_X_noint(rbind(mZ1, Z_i_demean), n_lags)
    }

    ## Produce forecasts


    ## Storage
    if (i %% thin == 0) {

      if (n_fcst > 0) {
        mu <- c(startpara$mu, numeric(n_fac))
        phi <- startpara$phi
        sigma <- startpara$sigma
        volatility_pred <- startlatent[TT, ]

        Z_pred <- matrix(0, n_fcst+n_lags, n_vars)
        Z_pred[1:n_lags, ] <- Z_i_demean[(TT-n_lags+1):TT,]
        for (j in 1:n_fcst) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- mfbvar:::create_X_t_noint(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
        Z_fcst[,,i/thin] <- Z_pred + rbind(d[(TT+1):(TT+n_lags), ,drop = FALSE], d_fcst[,,drop=FALSE]) %*% t(matrix(psi_i, nrow = n_vars))
      }

      Pi[,,i/thin] <- Pi_i
      Z[,,i/thin] <- Z_i
      psi[i/thin, ] <- psi_i

      mu_storage[,i/thin] <- startpara$mu
      sigma_storage[,i/thin] <- startpara$sigma
      phi_storage[,i/thin] <- startpara$phi

      fac_storage[,,i/thin] <- startfac
      facload_storage[,,i/thin] <- startfacload

      latent[,,i/thin] <- startlatent
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

    U <- mfbvar:::build_U_cpp(Pi = Pi_i, n_determ = n_determ, n_vars = n_vars, n_lags = n_lags)
    post_psi_Omega <- mfbvar:::posterior_psi_Omega_fsv(U = U, D_mat = D_mat, idivar = exp(startlatent[, 1:n_vars]),
                                          prior_psi_Omega = prior_psi_Omega)
    Y_tilde <- Z_i - tcrossprod(X, Pi_i) - t(startfac) %*% t(startfacload)

    post_psi <- mfbvar:::posterior_psi_mean_fsv(U = U, D_mat = D_mat, idivar = exp(startlatent[, 1:n_vars]), prior_psi_Omega = prior_psi_Omega,
                                   post_psi_Omega = post_psi_Omega, Y_tilde = Y_tilde, prior_psi_mean = prior_psi_mean)
    psi_i <- t(mfbvar:::rmultn(m = post_psi, Sigma = post_psi_Omega))

    mu_mat <- dt %*% t(matrix(psi_i, nrow = n_vars))

    if (verbose) {
      pb$tick()
    }
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Z = Z, psi = psi, Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst,
                     prior_Pi_Omega = prior_Pi_Omega, d = d, Y = x$Y, n_T = n_T, n_T_ = TT, n_reps = n_reps-1,
                     n_determ = n_determ, facload = facload_storage, latent = latent, mu = mu, sigma = sigma, phi = phi,
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
