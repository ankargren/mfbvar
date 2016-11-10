gibbs_sampler <- function(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, n_reps, n_fcst = NULL, lH, check_roots = FALSE,
                          init_Pi = NULL, init_Sigma = NULL, init_psi = NULL, init_Z = NULL,
                          d_fcst = NULL) {
  # prior_pi: (p * kp) matrix of prior mean for Pi
  # prior_pi_omega: (kp^2 * kp^2) matrix of prior covariance matrix for vec(pi')
  # prior_nu: scalar with the prior for nu
  # prior_s: (p * p) matrix with the prior for s
  # prior_psi:
  # prior_psi_omega:

  # Y: (T * p) matrix of main data (with NA where observations are missing)
  # d: (T * m) matrix of deterministic data
  # Lambda:


  # n_burnin: scalar with the number of burn-in replications
  # n_lags: scalar with the number of lags
  # n_reps: scalar with the number of replications

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_pi)))/n_vars^2
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_lags

  ################################################################
  # Preallocation
  # Pi and Sigma store their i-th draws in the third dimension, psi
  # is vectorized so it has its i-th draw stored in the i-th row
  # Pi:    p * pk * n_reps, each [,,i] stores Pi'
  # Sigma: p * p  * n_reps
  # psi:   n_reps * p
  # pre_Z: (T_ + k + 1) * p * n_reps
  # Z:     T * p * n_reps
  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps))
  psi   <- array(NA, dim = c(n_reps, n_vars * n_determ))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  if (!is.null(n_fcst)) {
    Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                   dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
    d_fcst_lags <- matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst), nrow = n_fcst + n_lags)
  }
  roots <- vector("numeric", n_reps)
  num_tries <- roots

  ################################################################
  # Gibbs sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the Gibbs sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  if (is.null(init_Z)) {
    Z[,, 1] <- fill_na(Y)
  } else {
    Z[,, 1] <- init_Z
  }

  ols_results <- ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ)

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    Pi[,, 1]    <- init_Pi
  }

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    Sigma[,, 1] <- init_Sigma
  }

  if (is.null(init_psi)) {
    psi[1, ] <- ols_results$psi
  } else {
    psi[1, ] <- init_psi
  }


  ################################################################
  # Create D (does not vary in the sampler), and find roots of Pi
  # if requested
  D <- build_DD(d = d, n_lags = n_lags)
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,,1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- max_eig_cpp(Pi_comp)
  }


  ################################################################
  # Calculations for the simulation smoother
  # These are the terms that do not vary in the sampler

  lH0 <- vector("list", n_T_)
  for(iter in 1:n_T_) {
    lH0[[iter]] = matrix(lH[[iter]][!is.na(Y[n_lags + iter,]),], ncol = n_vars * n_lags)
  }

  for (r in 2:(n_reps)) {

    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,, r-1] - d %*% t(matrix(psi[r-1, ], nrow = n_vars))
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)

    ################################################################
    ### Pi and Sigma step

    # Preliminary calculations
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    pi_sample <- solve(crossprod(XX)) %*% crossprod(XX, YY)

    # Posterior moments of Pi
    post_pi_omega <- solve(solve(prior_pi_omega) + crossprod(XX))
    post_pi       <- post_pi_omega %*% (solve(prior_pi_omega) %*% prior_pi + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% pi_sample)
    pi_diff <- prior_pi - pi_sample
    post_s <- prior_s + s_sample + t(pi_diff) %*% solve(post_pi_omega + solve(crossprod(XX)), pi_diff)
    nu <- n_T + prior_nu # Is this the right T? Or should it be T - lags?
    Sigma[,,r] <- rinvwish(v = nu, S = post_s)


    # Draw Pi conditional on Sigma
    # This ensures that the draw is stationary
    stationarity_check <- FALSE
    iter <- 0
    Pi_temp <- array(NA, dim = c(n_vars, n_vars * n_lags, 1000))
    while(stationarity_check == FALSE) {
      iter <- iter + 1
      Pi_temp[,,iter] <- rmatn(M = post_pi, Q = Sigma[,,r], P = post_pi_omega)
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
    # Steady-state step
    U <- build_U_cpp(Pi = Pi[,,r], n_determ = n_determ,
                     n_vars = n_vars, n_lags = n_lags)
    post_psi_omega <- posterior_psi_omega(U = U, D = D, sigma = Sigma[,, r],
                                          prior_psi_omega = prior_psi_omega)
    Y_tilde <- build_Y_tilde(Pi = Pi[,, r], z = Z[,, r-1])

    post_psi <- posterior_psi(U = U, D = D, sigma = Sigma[,, r], prior_psi_omega = prior_psi_omega,
                              psi_omega = post_psi_omega, Y_tilde = Y_tilde, prior_psi = prior_psi) # Seems to be correct
    psi[r, ] <- t(rmultn(m = post_psi, Sigma = post_psi_omega)) # Seems to be correct



    ################################################################
    # Smoothing step
    Q_comp     <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
    Q_comp[1:n_vars, 1:n_vars] <- t(chol(Sigma[,,r]))

    # Demean before putting into simulation smoother using the most recent draw of psi
    mZ <- Y - d %*% t(matrix(psi[r, ], nrow = n_vars))
    mZ <- mZ[-(1:n_lags), ]
    demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(psi[r-1, ], nrow = n_vars))
    h0 <- matrix(t(demeaned_z0), ncol = 1)
    h0 <- h0[(n_vars*n_lags):1,,drop = FALSE] # have to reverse the order

    smoothed_Z <- smooth_samp(mZ = mZ, mX = matrix(0, n_T_), lH = lH, lH0 = lH0,
                              mF = Pi_comp, mB = matrix(0, 16), mQ = Q_comp, iT = n_T_,
                              ip = n_lags, iq  = n_lags * n_vars, is = 1, h0 = h0, P0 = NULL,
                              X0 = 0)
    # For now, I'm just inserting h0 in the beginning of Z. Right now, Z has n_T number of rows,
    # but smooth_samp puts out n_T - n_lags rows.
    Z[,, r] <- rbind(demeaned_z0, smoothed_Z$mh[, 1:n_vars]) +
      d %*% t(matrix(psi[r, ], nrow = n_vars))


    ################################################################
    # Forecasting step


    if (!is.null(n_fcst)) {

      # Forecast the process with mean subtracted
      Z_fcst[1:n_lags, , r] <- Z[(n_T - n_lags+1):n_T,, r] - d[(n_T - n_lags+1):n_T, ] %*% t(matrix(psi[r, ], nrow = n_vars))
      for (h in 1:n_fcst) {
        Z_fcst[n_lags + h, , r] <- Pi[,, r] %*% matrix(c(t(Z_fcst[(n_lags+h-1):h,, r])), ncol = 1) +
          rmultn(m = matrix(0, nrow = 4), Sigma = Sigma[,,r])
      }

      # Add the mean
      Z_fcst[, , r] <- Z_fcst[, , r] + d_fcst_lags %*% t(matrix(psi[r, ], nrow = n_vars))
    }

  }

  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z)
  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (!is.null(n_fcst)) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)

}





