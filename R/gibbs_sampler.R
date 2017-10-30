#' Gibbs sampler for Mixed-Frequency BVAR
#'
#' \code{gibbs_sampler} runs a Gibbs sampler to approximate the posterior of the VAR model parameters.
#'
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar d_fcst TRUE
#' @templateVar Lambda TRUE
#' @templateVar prior_Pi_mean TRUE
#' @templateVar prior_Pi_Omega TRUE
#' @templateVar prior_S TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_reps TRUE
#' @templateVar init_Pi TRUE
#' @templateVar init_Sigma TRUE
#' @templateVar init_psi TRUE
#' @templateVar init_Z TRUE
#' @templateVar smooth_state TRUE
#' @templateVar check_roots TRUE
#' @templateVar verbose TRUE
#' @template man_template
#'
#' @details
#' The prior covariance of \eqn{\Pi} given \eqn{\Sigma} is \eqn{\Sigma \otimes \Omega_\Pi}.
#'
#' The function \code{gibbs_sampler_qf} is a temporary version of \code{gibbs_sampler} which is customized for quarterly data (for improved speed). Its purpose is to accomomdate conditional forecasting. The function \code{gibbs_sampler2} uses an R implementation of the simulation smoother, whereas \code{gibbs_sampler} uses a C++ implementation.
#'
#' @return
#' An object of class mfbvar.

gibbs_sampler <- function(Y, d, d_fcst = NULL, Lambda, prior_Pi_mean, prior_Pi_Omega, prior_S, prior_nu, prior_psi_mean, prior_psi_Omega,
                          n_fcst = NULL, n_reps, init_Pi = NULL, init_Sigma = NULL, init_psi = NULL, init_Z = NULL, smooth_state = FALSE, check_roots = TRUE,
                          verbose = TRUE) {

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_Pi_mean)))/n_vars^2
  n_pseudolags <- dim(Lambda)[2]/n_vars
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags

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
    Z_fcst_new <- Z_fcst
    d_fcst_lags <- matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst), nrow = n_fcst + n_lags)
  }
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  }
  if (smooth_state == TRUE) {
    smoothed_Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
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
    if (all(dim(Z[,, 1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,, 1]), collapse = " x ")))
    }

  }

  ols_results <- tryCatch(ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ),
                          error = function(cond) NULL)
  if (is.null(ols_results)) {
    ols_results <- list()
    ols_results$Pi <- prior_Pi_mean
    ols_results$S <- prior_S
    ols_results$psi <- prior_psi_mean
  }

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,, 1], n_vars = n_vars, n_lags = n_lags)
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
  D_mat <- build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  Z_1 <- Z[1:n_pseudolags,, 1]

  Pi_new <- Pi
  Sigma_new <- Sigma
  psi_new <- psi
  Z_new <- Z
  num_tries_new <- num_tries
  roots_new <- roots
  for (r in 2:(n_reps)) {
    ################################################################
    ### Pi and Sigma step
    #(Z_r1,             d,     psi_r1,                            prior_Pi_mean, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    set.seed(r)
    Pi_Sigma <- posterior_Pi_Sigma(Z_r1 = Z[,, r-1], d = d, psi_r1 = psi[r-1, , drop = FALSE], prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi[,,r]      <- Pi_Sigma$Pi_r
    Sigma[,,r]   <- Pi_Sigma$Sigma_r
    num_tries[r] <- Pi_Sigma$num_try
    roots[r]     <- Pi_Sigma$root

    set.seed(r)
    Pi_Sigma_new <- posterior_Pi_Sigma(Z_r1 = Z_new[,, r-1], d = d, psi_r1 = psi_new[r-1, , drop = FALSE], prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi_new[,,r]      <- Pi_Sigma_new$Pi_r
    Sigma_new[,,r]   <- Pi_Sigma_new$Sigma_r
    num_tries_new[r] <- Pi_Sigma_new$num_try
    roots_new[r]     <- Pi_Sigma_new$root


    ################################################################
    ### Steady-state step
    #(Pi_r,            Sigma_r,               Z_r1,             prior_psi_mean, prior_psi_Omega, D, n_vars, n_lags, n_determ)
    set.seed(r)
    psi[r, ] <- posterior_psi(Pi_r = Pi[,, r], Sigma_r = Sigma[,, r], Z_r1 = Z[,, r-1], prior_psi_mean, prior_psi_Omega, D_mat, n_vars, n_lags, n_determ)
    set.seed(r)
    psi_new[r, ] <- posterior_psi(Pi_r = Pi_new[,, r], Sigma_r = Sigma_new[,, r], Z_r1 = Z_new[,, r-1], prior_psi_mean, prior_psi_Omega, D_mat, n_vars, n_lags, n_determ)

    ################################################################
    ### Smoothing step
    #(Y, d, Pi_r,                                                             Sigma_r,               psi_r,                          Z_1, Lambda, n_vars, n_lags,                n_T_, smooth_state)
    set.seed(r)
    Z_res <- posterior_Z(Y, d, Pi_r = cbind(Pi[,, r], matrix(0, n_vars, n_vars*(n_pseudolags - n_lags))), Sigma_r = Sigma[,, r], psi_r = psi[r, , drop = FALSE], Z_1, Lambda, n_vars, n_lags = n_pseudolags, n_T_, smooth_state)
    Z[,, r] <- Z_res$Z_r

    Pi_r <- cbind(Pi[,,r], 0)
    mZ <- Y - d %*% t(matrix(psi[r,], nrow = n_vars))
    mZ <- as.matrix(mZ)
    demeaned_z0 <- Z_1 - d[1:n_lags, ] %*% t(matrix(psi[r,], nrow = n_vars))
    n_q <- 1
    T_b <- max(which(!apply(apply(Y[, 1:4, drop = FALSE], 2, is.na), 1, any)))
    Lambda_ <- build_Lambda(rep("average", n_q), 3)

    set.seed(r)
    Z_res_new <- kf_sim_smooth(mZ, cbind(Pi[,,r], 0), Sigma[,,r], Lambda_, demeaned_z0, n_q, T_b)
    Z_new_draw <- rbind(demeaned_z0, Z_res_new$Z_draw) + d %*% t(matrix(psi[r,], nrow = n_vars))
    Z_new[,, r] <- Z_new_draw

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

      Z_fcst_new[1:n_lags, , r] <- Z_new[(n_T - n_lags+1):n_T,, r] - d[(n_T - n_lags+1):n_T, ] %*% t(matrix(psi_new[r, ], nrow = n_vars))
      for (h in 1:n_fcst) {
        Z_fcst_new[n_lags + h, , r] <- Pi_new[,, r] %*% matrix(c(t(Z_fcst_new[(n_lags+h-1):h,, r])), ncol = 1) +
          rmultn(m = matrix(0, nrow = n_vars), Sigma = Sigma_new[,,r])
      }

      # Add the mean
      Z_fcst_new[, , r] <- Z_fcst_new[, , r] + d_fcst_lags %*% t(matrix(psi_new[r, ], nrow = n_vars))
    }

  }

  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, mdd = NULL, smoothed_Z = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = prior_nu, post_nu = n_T + prior_nu, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = prior_psi_Omega, prior_psi_mean = prior_psi_mean, n_reps = n_reps, Lambda = Lambda,
                     Pi_new = Pi_new, Sigma_new = Sigma_new, psi_new = psi_new, Z_new = Z_new, roots_new = roots_new, num_tries_new = num_tries_new,
                     Z_fcst_new = NULL)

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (!is.null(n_fcst)) {
    return_obj$Z_fcst <- Z_fcst
    return_obj$Z_fcst_new <- Z_fcst_new
  }
  if (smooth_state == TRUE) {
    return_obj$smoothed_Z <- smoothed_Z
  }

  return(return_obj)

}

