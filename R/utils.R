# The below loop just gets the error variances from AR(4) regressions

compute_error_variances <- function(Y) {
  n_vars <- ncol(Y)
  error_variance <- rep(NA, n_vars)
  for (i in 1:n_vars) {
    success <- NULL
    init_order <- 4
    for (ar_order in init_order:1) {
      error_variance[i] <- tryCatch(arima(na.omit(Y[,i]), order = c(ar_order, 0, 0), method = "ML")$sigma2,
                                    error = function(cond) NA)
      if (!is.na(error_variance[i])) {
        break
      } else {
        if (init_order < 1) {
          error_variance[i] <- var(na.omit(Y[,i]))
        }
      }
    }
  }
  return(error_variance)
}

check_required_params <- function(x, ...) {
  required_params <- unlist(list(...))
  test <- vapply(x[required_params], is.null, logical(1))
  if (any(test)) {
    stop("Missing elements: ", paste(required_params[which(test)], collapse = " "))
  }
}

list_to_variables <- function(x, envir, ...) {
  variables <- unlist(list(...))
  foo <- mapply(function(x1, x2, envir) assign(x1, x2, envir), names(x[variables]), x[variables], MoreArgs = list(envir = envir))
  invisible(NULL)
}

variable_initialization <- function(Y, freq, freqs, n_lags, Lambda_, n_thin) {
  n_vars <- ncol(Y)
  n_q <- sum(freq == freqs[1])
  if (n_q < n_vars) {
    T_b <- max(which(!apply(apply(Y[, freq == freqs[2], drop = FALSE], 2, is.na), 1, any)))
  } else {
    T_b <- nrow(Y)
  }

  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  n_thin <- ifelse(is.null(n_thin), 1, n_thin)

  Z_1 <- Y[1:n_pseudolags, ]
  return(list(n_vars = n_vars, n_q = n_q, T_b = T_b,
              n_pseudolags = n_pseudolags, n_T = n_T, n_T_ = n_T_,
              n_thin = n_thin, Z_1 = Z_1))
}

parameter_initialization <- function(Y, n_vars, n_lags, n_T_, init,
                                     n_fac = NULL, n_determ = NULL,
                                     n_sv = NULL,
                                     fsv, csv, ...) {
  arguments <- list(...)
  parameters <- unlist(arguments)
  steady_state <- "psi" %in% parameters
  if (fsv) {
    error_variance <- mfbvar:::compute_error_variances(Y)
  }

  const_latent <-  if (fsv) c(log(error_variance), rep(0, n_fac)) else 0

  init_available <- paste0("init_", parameters) %in% names(init)
  init_required <- parameters[!init_available]

  if (any(init_available)) {
    list_to_variables(init, parent.frame(), paste0("init_", parameters)[init_available])
  }

  for (i in seq_along(init_required)) {
    initval <- switch(init_required[i],
                      Z = mfbvar:::fill_na(Y),
                      psi = colMeans(mfbvar:::fill_na(Y)),
                      Pi = matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags)+!steady_state),
                      Sigma = cov(mfbvar:::fill_na(Y)),
                      omega = ifelse(!is.null(arguments$prior_psi_Omega),
                                     diag(prior_psi_Omega),
                                     rep(0.1, n_determ*n_vars)),
                      phi_mu = 1,
                      lambda_mu = 1,
                      mu = log(error_variance),
                      sigma = rep(0.75, n_sv),
                      phi = rep(0.2, n_sv),
                      facload = matrix(rnorm(n_vars*n_fac, sd = 0.5)^2, n_vars, n_fac),
                      f = matrix(rnorm(n_fac * n_T_, sd = 0.5), n_fac, n_T_),
                      latent = t(cbind(matrix(const_latent, nrow = n_T_, ncol = n_sv, byrow = TRUE))),
                      latent0 = numeric(n_vars + n_fac),
                      global = 0.1,
                      aux = rep(0.1, n_vars^2 * n_lags),
                      local = rep(0.1, n_vars^2 * n_lags),
                      slice =  rep(1, n_vars^2 * n_lags)
    )
    assign(paste0("init_", init_required[i]), initval)
  }

  return(mget(paste0("init_", parameters)))
}

storage_initialization <- function(init_params, params, envir, n_vars, n_lags,
                                  n_reps, n_thin, n_T, n_T_, n_determ = NULL,
                                  n_fac = NULL, n_fcst, n_sv = NULL) {
  steady_state <- "psi" %in% params

  for (i in seq_along(params)) {
    initval <- init_params[[paste0("init_", params[i])]]
    assign(params[i],
           switch(params[i],
    Pi = array(initval, dim = c(n_vars, n_vars*n_lags+!steady_state, n_reps/n_thin)),
    Sigma = array(initval, dim = c(n_vars, n_vars, n_reps/n_thin)),
    psi = array(initval, dim = c(n_vars * n_determ, 1, n_reps/n_thin)),
    Z = array(initval, dim = c(n_T, n_vars, n_reps/n_thin)),
    mu = array(initval, dim = c(n_vars, 1, n_reps/n_thin)),
    sigma = array(initval, dim = c(n_sv, 1, n_reps/n_thin)),
    phi = array(initval, dim = c(n_sv, 1, n_reps/n_thin)),
    facload = array(matrix(initval, nrow = n_vars, ncol = n_fac),
                     dim = c(n_vars, n_fac, n_reps/n_thin)),
    f = array(matrix(initval, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin)),
    latent = array(initval, dim = c(n_T_, n_sv, n_reps/n_thin),
                    dimnames = list(rownames(initval), colnames(initval), NULL)),
    omega = array(initval, dim = c(n_vars * n_determ, 1, n_reps/n_thin)),
    phi_mu = array(initval, dim = c(1, 1, n_reps/n_thin)),
    lambda_mu = array(initval, dim = c(1, 1, n_reps/n_thin)),
    aux = array(initval, dim = c(n_vars*n_vars*n_lags, 1, n_reps/n_thin)),
    local = array(initval, dim = c(n_vars*n_vars*n_lags, 1, n_reps/n_thin)),
    global = array(initval, dim = c(1, 1, n_reps/n_thin)),
    slice = array(initval, dim = c(n_vars*n_vars*n_lags, 1, n_reps/n_thin))),
    envir)
  }

  Z_fcst <- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps/n_thin))
  if (n_fcst > 0) {
    rownames(Z_fcst) <- c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst))
    Z_fcst[,,1] <- 0
  } else {
    rownames(Z_fcst) <- (n_T-n_lags+1):n_T
  }
  assign("Z_fcst", Z_fcst, envir)
}

fsv_initialization <- function(priorsigmaidi, priorsigmafac, priormu,
                               priorfacload, restrict, priorphiidi, priorphifac,
                               n_vars, n_fac) {

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
  return(list(priorsigmaidi = priorsigmaidi,
              priorsigmafac = priorsigmafac,
              bmu = bmu,
              Bmu = Bmu,
              Bsigma = Bsigma,
              B011inv = B011inv,
              B022inv = B022inv,
              armatau2 = armatau2,
              armarestr = armarestr,
              a0idi = a0idi,
              b0idi = b0idi,
              a0fac = a0fac,
              b0fac = b0fac,
              priorh0 = priorh0,
              n_sv = n_fac + n_vars))
}

ssng_initialization <- function(prior_ng, s) {
  c0 <- ifelse(is.null(prior_ng), 0.01, prior_ng[1])
  c1 <- ifelse(is.null(prior_ng), 0.01, prior_ng[2])
  s <- ifelse(is.null(s), 1, s)
  return(list(c0 = c0, c1 = c1, s = s))
}

ss_initialization <- function(d, d_fcst, n_T, n_lags, n_fcst) {
  n_determ <- if (!is.null(d)) dim(d)[2] else NULL
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]
  return(list(d_fcst_lags = d_fcst_lags, D_mat = D_mat, dt = dt, d1 = d1,
              n_determ = n_determ))
}

dl_initialization <- function(a, gig, n_cores) {
  a   <- ifelse(is.null(a), 1, a)
  gig <- ifelse(is.null(gig), TRUE, FALSE)
  RcppParallel::setThreadOptions(numThreads = x$n_cores)
}


csv_initialization <- function(prior_phi, prior_sigma2) {
  phi_invvar <- 1/prior_phi[2]
  phi_meaninvvar <- prior_phi[1] * phi_invvar
  prior_df <- prior_sigma2[2]
  prior_sigma2 <- prior_sigma2[1]

  return(list(phi_invvar = phi_invvar, phi_meaninvvar = phi_meaninvvar,
              prior_sigma2 = prior_sigma2, prior_df = prior_df,
              n_sv = 1))
}
