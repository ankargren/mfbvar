#' @description Older, slower version of the Gibbs sampler.
#' @inherit gibbs_sampler
#' @templateVar lH TRUE
#' @template man_template
#'
gibbs_sampler2 <- function(prior_Pi_mean, prior_Pi_Omega, prior_nu, prior_S, prior_psi_mean, prior_psi_Omega,
                           Y, d, n_reps, n_fcst = NULL, lH, check_roots = TRUE,
                           init_Pi = NULL, init_Sigma = NULL, init_psi = NULL, init_Z = NULL,
                           d_fcst = NULL, smooth_state = FALSE) {

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_Pi_mean)))/n_vars^2
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
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  if (!is.null(n_fcst)) {
    Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                   dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
    d_fcst_lags <- matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst), nrow = n_fcst + n_lags)
  }
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  }
  if (smooth_state == TRUE) {
    smoothed_Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
    smoothed_Y     <- array(NA, dim = c(n_T, n_vars, n_reps))
  }



  ################################################################
  ### Gibbs sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the Gibbs sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  if (is.null(init_Z)) {
    Z[,, 1] <- fill_na(Y)
  } else {
    if (all(dim(Z[,,1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,,1]), collapse = " x ")))
    }

  }

  ols_results <- ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ)

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,,1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,,1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,,1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- max_eig_cpp(Pi_comp)
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
    if (roots[1] < 1) {
      psi[1, ] <- ols_results$psi
    } else {
      psi[1, ] <- prior_psi_mean
    }
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

  # For the posterior of Pi
  Omega_Pi <- solve(prior_Pi_Omega) %*% prior_Pi_mean

  # Calculations for the simulation smoother
  lH0 <- vector("list", n_T_)
  for(iter in 1:n_T_) {
    lH0[[iter]] = matrix(lH[[iter]][!is.na(Y[n_lags + iter,]),], ncol = n_vars * n_lags)
  }
  mX <- matrix(0, n_T_)
  mB <- matrix(0, n_vars*n_lags)


  for (r in 2:(n_reps)) {
    ################################################################
    ### Preliminary calculations

    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,, r-1] - d %*% t(matrix(psi[r-1, ], nrow = n_vars))
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    Pi_sample <- solve(crossprod(XX)) %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_Pi_Omega <- solve(solve(prior_Pi_Omega) + crossprod(XX))
    post_Pi       <- post_Pi_Omega %*% (Omega_Pi + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% Pi_sample)
    Pi_diff <- prior_Pi_mean - Pi_sample
    post_s <- prior_S + s_sample + t(Pi_diff) %*% solve(prior_Pi_Omega + solve(crossprod(XX))) %*% Pi_diff
    nu <- n_T + prior_nu # Is this the right T? Or should it be T - lags?
    Sigma[,,r] <- rinvwish(v = nu, S = post_s)


    # Draw Pi conditional on Sigma
    # This ensures that the draw is stationary
    stationarity_check <- FALSE
    iter <- 0
    Pi_temp <- array(NA, dim = c(n_vars, n_vars * n_lags, ifelse(check_roots, 1000, 1)))
    while(stationarity_check == FALSE) {
      iter <- iter + 1
      Pi_temp[,,iter] <- rmatn(M = t(post_Pi), Q = post_Pi_Omega, P = Sigma[,,r])
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
    post_psi_Omega <- posterior_psi_Omega(U = U, D_mat = D, Sigma = Sigma[,, r],
                                          prior_psi_Omega = prior_psi_Omega)
    Y_tilde <- build_Y_tilde(Pi = Pi[,, r], z = Z[,, r-1])

    post_psi <- posterior_psi(U = U, D_mat = D, Sigma = Sigma[,, r], prior_psi_Omega = prior_psi_Omega,
                              post_psi_Omega = post_psi_Omega, Y_tilde = Y_tilde, prior_psi_mean = prior_psi_mean)
    psi[r, ] <- t(rmultn(m = post_psi, Sigma = post_psi_Omega))



    ################################################################
    ### Smoothing step
    Q_comp     <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
    Q_comp[1:n_vars, 1:n_vars] <- t(chol(Sigma[,,r]))

    # Demean before putting into simulation smoother using the most recent draw of psi
    mZ <- Y - d %*% t(matrix(psi[r, ], nrow = n_vars))
    mZ <- as.matrix(mZ[-(1:n_lags), ])
    demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(psi[r, ], nrow = n_vars))
    h0 <- matrix(t(demeaned_z0[n_lags:1,]), ncol = 1)

    simulated_Z <- smooth_samp_u3(mZ = mZ, mX = mX, lH = lH, lH0 = lH0,
                                  mF = Pi_comp, mB = mB, mQ = Q_comp, iT = n_T_,
                                  ip = n_vars, iq  = n_lags * n_vars, is = 1, h0 = h0, P0 = NULL,
                                  X0 = 0)
    # For now, I'm just inserting h0 in the beginning of Z. Right now, Z has n_T number of rows,
    # but smooth_samp puts out n_T - n_lags rows.
    Z[,, r] <- rbind(demeaned_z0, simulated_Z$mh[, 1:n_vars]) +
      d %*% t(matrix(psi[r, ], nrow = n_vars))

    # Also save the smoothed value of the state
    if (smooth_state == TRUE) {
      smoothed_state <- smoothing_state(mZ = mZ, mX = mX, lH = lH, lH0 = lH0,
                                        mF = Pi_comp, mB = mB, mQ = Q_comp, iT = n_T_,
                                        ip = n_vars, iq  = n_lags * n_vars, is = 1, h0 = h0, P0 = NULL,
                                        X0 = 0)
      smoothed_Z[,, r] <- rbind(demeaned_z0, smoothed_state$mh[, 1:n_vars]) +
        d %*% t(matrix(psi[r, ], nrow = n_vars))
      smoothed_Y[,, r] <- rbind(demeaned_z0, smoothed_state$mZ[, 1:n_vars]) +
        d %*% t(matrix(psi[r, ], nrow = n_vars))
    }

    ################################################################
    ### Forecasting step
    if (!is.null(n_fcst)) {

      # Forecast the process with mean subtracted
      Z_fcst[1:n_lags, , r] <- Z[(n_T - n_lags+1):n_T,, r] - d[(n_T - n_lags+1):n_T, ] %*% t(matrix(psi[r, ], nrow = n_vars))
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
                     Z_fcst = NULL, mdd = NULL, smoothed_Z = NULL, smoothed_Y = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = prior_nu, nu = nu, d = d, Y = Y, n_T = n_T, n_T_ = n_T_, lH0 = lH0,
                     prior_psi_Omega = prior_psi_Omega, prior_psi_mean = prior_psi_mean, n_reps = n_reps)

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (!is.null(n_fcst)) {
    return_obj$Z_fcst <- Z_fcst
  }
  if (smooth_state == TRUE) {
    return_obj$smoothed_Z <- smoothed_Z
    return_obj$smoothed_Y <- smoothed_Y
  }

  return(return_obj)

}
