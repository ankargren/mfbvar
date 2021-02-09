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

variable_initialization <- function(Y, freq, freqs, n_lags, Lambda_, d = NULL, d_fcst = NULL) {
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
    Y <- Y[complete_quarters, ]
    if (!is.null(d)) {
      d_fcst <- rbind(d[!complete_quarters, , drop = FALSE], d_fcst)
      d <- d[complete_quarters, , drop = FALSE]
    }
  }

  n_pseudolags <- max(c(n_lags, ncol(Lambda_)/nrow(Lambda_)))
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

  n_thin <- ifelse(is.null(n_thin), 1, n_thin)
  return(list(n_vars = n_vars, n_determ = n_determ, n_q = n_q, T_b = T_b,
              n_pseudolags = n_pseudolags, n_T = n_T, n_T_ = n_T_,
              n_thin = n_thin))
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
    assign(paste0("init_", init_required[i]), initval, parent.frame())
  }

  return(mget(paste0("init_", parameters)))
}
