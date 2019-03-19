mcmc_sampler.mfbvar_minn_sv <- function(x, ...){

  add_args <- list(...)

  n_reps <- add_args$n_reps
  thin <- x$thin

  init <- add_args$init

  init_Pi <- init$init_Pi
  init_Z <- init$init_Z

  init_mu <- init$init_mu
  init_sigma <- init$init_sigma
  init_phi <- init$init_phi

  init_facload <- init$facload
  init_fac <- init$fac

  init_latent <- init$init_latent
  init_latent0 <- init$init_latent0

  cl <- x$cl
  Z1 <- x$Z1

  n_vars <- nrow(init_Pi)
  n_lags <- (ncol(init_Pi)-1)/n_vars
  n_q <- x$n_q
  n_m <- n_vars - n_q
  n_fac <- x$n_fac
  n_fcst <- x$n_fcst

  y_in_p <- x$Y[-(1:n_lags), ]

  T_b <- min(apply(y_in_p[,1:n_m], 2, function(x) ifelse(any(is.na(x)), min(which(is.na(x))), Inf))-1, nrow(y_in_p))
  TT <- nrow(x$Y) - n_lags
  n_T <- nrow(x$Y)

  mf <- TRUE

  ## Set up cluster (if used)
  if (!is.null(cl)) {
    parallelize <- TRUE
    parallel::clusterCall(cl, fun = function() library(mfbvar))
    parallel::clusterExport(cl, varlist = c("par_fun"))
  } else {
    parallelize <- FALSE
  }

  prior_zero_mean <- all(x$prior_Pi_AR1) == 0
  if (prior_zero_mean) {
    if (n_vars*n_lags > 1.05 * TT) {
      par_fun <- par_fun_top(rmvn_bcm)
    } else {
      par_fun <- par_fun_top(rmvn_rue)
    }
  }

  ## Obtain the aggregation matrix for the quarterly only
  if (mf) {
    Lambda_companion <- build_Lambda(c(rep("m", n_vars-n_q), rep("q", n_q)), n_lags)
    Lambda_comp <- matrix(Lambda_companion[(n_m+1):n_vars, c(t(sapply((n_m+1):n_vars, function(x) seq(from = x, to = n_vars*n_lags, by = n_vars))))],
                          nrow = n_q)
  }

  Pi <- array(init_Pi, dim = c(n_vars, n_vars*n_lags + 1, n_reps/thin))
  Z <- array(y_in_p, dim = c(TT, n_vars, n_reps/thin))
  if (n_fcst > 0) {
    Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
            dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }

  mu_storage <- matrix(init_mu, n_vars, n_reps/thin)
  sigma_storage <- matrix(init_sigma, n_vars+n_fac, n_reps/thin)
  phi_storage <- matrix(init_phi, n_vars+n_fac, n_reps/thin)

  facload_storage <- array(matrix(init_facload, nrow = n_vars, ncol = n_fac), dim = c(n_vars, n_fac, n_reps/thin))
  fac_storage <- array(matrix(init_fac, n_fac, TT), n_reps/thin)

  latent <- array(init_latent, dim = c(TT, n_vars+n_fac, n_reps/thin),
                  dimnames = list(rownames(init_latent), colnames(init_latent), NULL))

  Pi_i <- Pi[,,1]
  Z_i <- init_Z
  startpara <- list(mu = init_mu,
                    phi = init_phi,
                    sigma = init_sigma)
  startlatent <- latent[,,1]
  startlatent0 <- init_latent0 #cbind(fLambda$idivol0, fLambda$facvol0) #
  startfacload <- init_facload #fLambda$facload #
  startfac <- init_fac #fLambda$f[,-(1:n_lags),drop=F] #

  prior_Pi_Omega <- x$prior_Pi_Omega

  error <- NULL
  for (i in 1:n_reps) {
    ## Square root of idiosyncratic variances (in dense form)
    Sig <- exp(0.5 * startlatent[, 1:n_vars])

    ## Mixed-frequency block: sample latent monthly series
    if (mf) {
      Z_i <- tryCatch(simsm_adaptive_univariate(y_in_p, Pi_i, Sig, Lambda_comp, Z1, n_q, T_b, t(startfac) %*% t(startfacload)), error = function(cond) cond)
      if (inherits(Z_i, "error")) {
        warning("MCMC halted because of an error in the mixed-frequency step. See $error for more information.")
        error <- list(error = Z_i, iter = i, block = "z")
        break
      }
      Z_i <- rbind(Z1, Z_i)
      X <- create_X(Z_i, n_lags)
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
        for (j in 1:(n_fcst+n_lags)) {
          volatility_pred <- mu + phi * (volatility_pred - mu) + rnorm(n_vars+n_fac, sd = sigma)
          error_pred <- rnorm(n_vars+n_fac, sd = exp(volatility_pred * 0.5))
          X_t <- create_X_t(Z_pred[j:(n_lags+j-1), ])
          Z_pred[j+n_lags, ] <- Pi_i %*% X_t + startfacload %*% error_pred[(n_vars+1):(n_vars+n_fac)] + error_pred[1:n_vars]
        }
      }

      Pi[,,i/thin] <- Pi_i
      Z[,,i/thin] <- Z_i
      Z_fcst[,,i/thin] <- Z_pred

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
                                                  startfacload = startfacload, startfac = startfac, ...), error = function(cond) cond)
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
          Pi_i[j,] <- tryCatch(rmvn(X/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5)), error = function(cond) cond)
        }
      } else {
        for (j in 1:n_vars) {
          Pi_i[j,] <- tryCatch(mvn_ccm(X/exp(startlatent[,j]*0.5), prior_Pi_Omega[,j], latent_nofac[,j]/exp(startlatent[,j]*0.5), prior_Pi_AR1[j], j), error = function(cond) cond)
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
  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = NULL, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = prior_nu, post_nu = prior_nu + n_T_, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = NULL, prior_psi_mean = NULL, n_reps = n_reps-1, Lambda = Lambda, freq = freq,
                     init = list(init_Pi = Pi[,, n_reps], init_Sigma = Sigma[,, n_reps], init_Z = Z[,, n_reps]))

  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}
