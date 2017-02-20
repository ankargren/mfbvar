#' @description The Gibbs sampler for quarterly data.
#' @inherit gibbs_sampler2
gibbs_sampler_qf <- function(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, n_reps, n_fcst = NULL, lH, check_roots = TRUE,
                          init_Pi = NULL, init_Sigma = NULL, init_psi = NULL, init_Z = NULL,
                          d_fcst = NULL, smooth_state = FALSE) {

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_pi)))/n_vars^2
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_lags

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

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps))
  psi   <- array(NA, dim = c(n_reps, n_vars * n_determ))
  Z     <- na.omit(Y)
  n_Z <- nrow(Z)
  d_full <- d
  non_NA <- which(apply(Y, 1, function(x) all(!is.na(x))))
  d <- d[non_NA,,drop=FALSE]
  last_complete <- max(non_NA)
  if (!is.null(n_fcst)) {
    Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                   dimnames = list(c((last_complete-n_lags+1):last_complete, paste0("fcst_", 1:n_fcst)), NULL, NULL))
    d_fcst_lags <- matrix(rbind(d_full[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst), nrow = n_fcst + n_lags)
  }
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  }


  ################################################################
  ### Gibbs sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the Gibbs sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  ols_results <- ols_initialization(z = Z, d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ)

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,,1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,,1]), collapse = " x ")))
    }
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
    psi[1, ] <- ols_results$psi
  } else {
    if (length(psi[1, ]) == length(init_psi)) {
      psi[1,] <- init_psi
    } else {
      stop(paste0("The length of init_psi is ", paste(length(init_psi), collapse = " x "), ", but should be ", paste(length(psi[1,]), collapse = " x ")))
    }
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  # if requested
  D <- build_DD(d = d, n_lags = n_lags)
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,,1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- max_eig_cpp(Pi_comp)
  }

  # For the posterior of Pi
  omega_pi <- solve(prior_pi_omega) %*% prior_pi

  for (r in 2:(n_reps)) {
    ################################################################
    ### Preliminary calculations

    # Demean z, create Z (companion form version)
    demeaned_z <- Z - d %*% t(matrix(psi[r-1, ], nrow = n_vars))
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    pi_sample <- solve(crossprod(XX)) %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_pi_omega <- solve(solve(prior_pi_omega) + crossprod(XX))
    post_pi       <- post_pi_omega %*% (omega_pi + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% pi_sample)
    pi_diff <- prior_pi - pi_sample
    post_s <- prior_s + s_sample + t(pi_diff) %*% solve(post_pi_omega + solve(crossprod(XX))) %*% pi_diff
    nu <- n_T + prior_nu # Is this the right T? Or should it be T - lags?
    Sigma[,,r] <- rinvwish(v = nu, S = post_s)


    # Draw Pi conditional on Sigma
    # This ensures that the draw is stationary
    stationarity_check <- FALSE
    iter <- 0
    Pi_temp <- array(NA, dim = c(n_vars, n_vars * n_lags, ifelse(check_roots, 1000, 1)))
    while(stationarity_check == FALSE) {
      iter <- iter + 1
      Pi_temp[,,iter] <- rmatn(M = t(post_pi), Q = post_pi_omega, P = Sigma[,,r])
      Pi_comp    <- build_companion(Pi_temp[,, iter], n_vars = n_vars, n_lags = n_lags)
      if (check_roots == TRUE) {
        roots[r] <- max_eig_cpp(Pi_comp)
      }
      if (roots[r] < 1) {
        stationarity_check <- TRUE
        num_tries[r] <- iter
        Pi[,,r] <- Pi_temp[,,iter]
      }
      if (iter == 1000) {
        stop("Attempted to draw stationary Pi 1,000 times.")
      }
    }


    ################################################################
    ### Steady-state step
    U <- build_U_cpp(Pi = Pi[,,r], n_determ = n_determ,
                     n_vars = n_vars, n_lags = n_lags)
    post_psi_omega <- posterior_psi_omega(U = U, D_mat = D, Sigma = Sigma[,, r],
                                          prior_psi_omega = prior_psi_omega)
    Y_tilde <- build_Y_tilde(Pi = Pi[,, r], z = Z)

    post_psi <- posterior_psi(U = U, D_mat = D, Sigma = Sigma[,, r], prior_psi_omega = prior_psi_omega,
                              post_psi_omega = post_psi_omega, Y_tilde = Y_tilde, prior_psi = prior_psi)
    psi[r, ] <- t(rmultn(m = post_psi, Sigma = post_psi_omega))


        ################################################################
    ### Forecasting step
    if (!is.null(n_fcst)) {

      # Forecast the process with mean subtracted
      Z_fcst[1:n_lags, , r] <- Z[(n_Z - n_lags+1):n_Z,] - d_full[(n_Z - n_lags+1):n_Z, ] %*% t(matrix(psi[r, ], nrow = n_vars))
      for (h in 1:n_fcst) {
        Z_fcst[n_lags + h, , r] <- Pi[,, r] %*% matrix(c(t(Z_fcst[(n_lags+h-1):h,, r])), ncol = 1) +
          rmultn(m = matrix(0, nrow = n_vars), Sigma = Sigma[,,r])
      }

      # Add the mean
      Z_fcst[, , r] <- Z_fcst[, , r] + d_fcst_lags %*% t(matrix(psi[r, ], nrow = n_vars))
    }

  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, mdd = NULL)

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (!is.null(n_fcst)) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}





