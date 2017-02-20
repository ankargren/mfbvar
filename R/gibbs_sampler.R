#' Gibbs sampler for Mixed-Frequency BVAR
#'
#' \code{gibbs_sampler} runs a Gibbs sampler to approximate the posterior of the VAR model parameters.
#'
#' @templateVar prior_pi TRUE
#' @templateVar prior_pi_omega TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_s TRUE
#' @templateVar prior_psi TRUE
#' @templateVar prior_psi_omega TRUE
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar n_reps TRUE
#' @templateVar n_fcst TRUE
#' @templateVar Lambda TRUE
#' @templateVar check_roots TRUE
#' @templateVar init_Pi TRUE
#' @templateVar init_Sigma TRUE
#' @templateVar init_psi TRUE
#' @templateVar init_Z TRUE
#' @templateVar d_fcst TRUE
#' @templateVar smooth_state TRUE
#' @template man_template
#'
#' @details
#' The prior covariance of \eqn{\Pi} given \eqn{\Sigma} is \eqn{\Sigma \otimes \Omega_\Pi}.
#'
#' The function \code{gibbs_sampler_qf} is a temporary version of \code{gibbs_sampler} which is customized for quarterly data (for improved speed). Its purpose is to accomomdate conditional forecasting. The function \code{gibbs_sampler2} uses an R implementation of the simulation smoother, whereas \code{gibbs_sampler} uses a C++ implementation.
#'
#' @return
#' An object of class mfbvar.

gibbs_sampler <- function(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, n_reps, n_fcst = NULL, Lambda, check_roots = TRUE,
                          init_Pi = NULL, init_Sigma = NULL, init_psi = NULL, init_Z = NULL,
                          d_fcst = NULL, smooth_state = FALSE) {

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_pi)))/n_vars^2
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
      psi[1, ] <- prior_psi
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
  inv_prior_pi_omega <- solve(prior_pi_omega)
  omega_pi <- inv_prior_pi_omega %*% prior_pi

  Z_1 <- Z[1:n_pseudolags,, 1]

  for (r in 2:(n_reps)) {
    ################################################################
    ### Pi and Sigma step
    #(Z_r1,             d,     psi_r1,                            prior_pi, inv_prior_pi_omega, omega_pi, prior_s, prior_nu, check_roots, n_vars, n_lags, n_T)
    pi_sigma <- pi_sigma_posterior(Z_r1 = Z[,, r-1], d = d, psi_r1 = psi[r-1, , drop = FALSE], prior_pi, inv_prior_pi_omega, omega_pi, prior_s, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi[,,r]      <- pi_sigma$Pi_r
    Sigma[,,r]   <- pi_sigma$Sigma_r
    num_tries[r] <- pi_sigma$num_try
    roots[r]     <- pi_sigma$root

    ################################################################
    ### Steady-state step
    #(Pi_r,            Sigma_r,               Z_r1,             prior_psi, prior_psi_omega, D, n_vars, n_lags, n_determ)
    psi[r, ] <- psi_posterior(Pi_r = Pi[,, r], Sigma_r = Sigma[,, r], Z_r1 = Z[,, r-1], prior_psi, prior_psi_omega, D_mat, n_vars, n_lags, n_determ)

    ################################################################
    ### Smoothing step
    #(Y, d, Pi_r,                                                             Sigma_r,               psi_r,                          Z_1, Lambda, n_vars, n_lags,                n_T_, smooth_state)
    Z_res <- Z_posterior(Y, d, Pi_r = cbind(Pi[,, r], matrix(0, n_vars, n_vars*(n_pseudolags - n_lags))), Sigma_r = Sigma[,, r], psi_r = psi[r, , drop = FALSE], Z_1, Lambda, n_vars, n_lags = n_pseudolags, n_T_, smooth_state)
    Z[,, r] <- Z_res$Z_r
    if (smooth_state == TRUE) {
      smoothed_Z[,, r] <- Z_res$smoothed_Z_r
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
                     Z_fcst = NULL, mdd = NULL, smoothed_Z = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, prior_pi_omega = prior_pi_omega, prior_pi = prior_pi,
                     prior_s = prior_s, prior_nu = prior_nu, nu = n_T + prior_nu, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_omega = prior_psi_omega, prior_psi = prior_psi, n_reps = n_reps, Lambda = Lambda)

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (!is.null(n_fcst)) {
    return_obj$Z_fcst <- Z_fcst
  }
  if (smooth_state == TRUE) {
    return_obj$smoothed_Z <- smoothed_Z
  }

  return(return_obj)

}

