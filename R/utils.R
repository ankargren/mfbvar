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

variable_initialization <- function(Y, freq, freqs, n_lags, Lambda_, n_thin,
                                    d = NULL, d_fcst = NULL) {
  n_vars <- ncol(Y)
  n_determ <- if (!is.null(d)) dim(d)[2] else NULL
  n_q <- sum(freq == freqs[1])
  if (n_q < n_vars) {
    T_b <- max(which(!apply(apply(Y[, freq == freqs[2], drop = FALSE], 2, is.na), 1, any)))
  } else {
    T_b <- nrow(Y)
  }
  if (n_q == 0 || n_q == n_vars) {
    complete_quarters <- apply(Y, 1, function(x) !any(is.na(x)))
    #Y <- Y[complete_quarters, ]
    if (!is.null(d)) {
      d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
      #d <- d[complete_quarters, , drop = FALSE]
    }
  }

  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  n_thin <- ifelse(is.null(n_thin), 1, n_thin)

  Z_1 <- Y[1:n_pseudolags, ]
  return(list(n_vars = n_vars, n_determ = n_determ, n_q = n_q, T_b = T_b,
              n_pseudolags = n_pseudolags, n_T = n_T, n_T_ = n_T_,
              n_thin = n_thin, Z_1 = Z_1))
}

parameter_initialization <- function(Y, n_vars, n_lags, n_T_, init,
                                     n_fac = NULL, n_determ = NULL, ...) {
  arguments <- list(...)
  parameters <- unlist(arguments)
  steady_state <- "psi" %in% parameters
  fsv <- "facload" %in% parameters
  if (fsv) {
    error_variance <- compute_error_variances(Y)
  }

  init_available <- paste0("init", parameters) %in% names(init)
  init_required <- parameters[!init_available]

  if (any(init_available)) {
    list_to_variables(init, parent.frame(), paste0("init", parameters)[init_available])
  }

  for (i in seq_along(init_required)) {
    initval <- switch(init_required[i],
                      Z = fill_na(Y),
                      psi = colMeans(fill_na(Y)),
                      Pi = matrix(0, nrow = n_vars, ncol = n_vars*(n_vars*n_lags)+!steady_state),
                      omega = ifelse(!is.null(arguments$prior_psi_Omega),
                                     diag(prior_psi_Omega),
                                     rep(0.1, n_determ*n_vars)),
                      phi_mu = 1,
                      lambda_mu = 1,
                      mu = log(error_variance),
                      sigma = rep(0.75, n_vars + n_fac),
                      phi = rep(0.2, n_vars + n_fac),
                      facload = matrix(rnorm(n_vars*n_fac, sd = 0.5)^2, n_vars, n_fac),
                      f = matrix(rnorm(n_fac * n_T_, sd = 0.5), n_fac, n_T_),
                      latent = t(cbind(matrix(c(log(error_variance), rep(1, n_fac)), nrow = n_T_, ncol = n_vars+n_fac, byrow = TRUE))),
                      latent0 <- numeric(n_vars + n_fac)
    )
    assign(paste0("init_", init_required[i]), initval)
  }

  return(mget(paste0("init_", parameters)))
}

storage_initialization <- function(init_params, params, envir, n_vars, n_lags,
                                  n_reps, n_thin, n_T, n_T_, n_determ = NULL,
                                  n_fac = NULL, n_fcst) {
  steady_state <- "psi" %in% params

  for (i in seq_along(params)) {
    initval <- init_params[[paste0("init_", params[i])]]
    assign(params[i],
           switch(params[i],
    Pi = array(initval, dim = c(n_vars, n_vars*n_lags+!steady_state, n_reps/n_thin)),
    psi = array(initval, dim = c(n_reps/n_thin, n_vars * n_determ)),
    Z = array(initval, dim = c(n_T, n_vars, n_reps/n_thin)),
    mu = matrix(initval, n_vars, n_reps/n_thin),
    sigma = matrix(initval, n_vars+n_fac, n_reps/n_thin),
    phi = matrix(initval, n_vars+n_fac, n_reps/n_thin),
    facload = array(matrix(initval, nrow = n_vars, ncol = n_fac),
                     dim = c(n_vars, n_fac, n_reps/n_thin)),
    f = array(matrix(initval, n_fac, n_T_), dim = c(n_fac, n_T_, n_reps/n_thin)),
    h = array(t(initval), dim = c(n_T_, n_vars+n_fac, n_reps/n_thin),
               dimnames = list(rownames(initval), colnames(initval), NULL)),
    omega = matrix(initval, nrow = n_reps/n_thin, ncol = n_vars * n_determ, byrow = TRUE),
    phi_mu = rep(initval, n_reps/n_thin),
    lambda_mu = rep(initval, n_reps/n_thin)),
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
              priorh0 = priorh0))
}

ssng_initialization <- function(prior_ng, s) {
  c0 <- ifelse(is.null(prior_ng), 0.01, prior_ng[1])
  c1 <- ifelse(is.null(prior_ng), 0.01, prior_ng[2])
  s <- ifelse(is.null(s), 1, s)
  return(list(c0 = c0, c1 = c1, s = s))
}

ss_initialization <- function(d, d_fcst, n_T, n_lags, n_fcst) {
  d_fcst_lags <- as.matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst))
  d_fcst_lags <- d_fcst_lags[1:(n_lags+n_fcst), , drop = FALSE]
  D_mat <- mfbvar:::build_DD(d = d, n_lags = n_lags)
  dt <- d[-(1:n_lags), , drop = FALSE]
  d1 <- d[1:n_lags, , drop = FALSE]
  return(list(d_fcst_lags = d_fcst_lags, D_mat = D_mat, dt = dt, d1 = d1))
}
