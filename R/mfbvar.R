mfbvar <- function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_mean_pi,
                   lambda1, lambda2, prior_nu = NULL, prior_mean_psi, prior_var_psi, ...) {

  fun_call <- match.call()
  n_vars <- ncol(Y)
  n_T <- nrow(Y)
  n_T_ <- n_T - n_lags

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_pi_sigma(lambda1 = lambda1, lambda2 = lambda2, prior_mean = prior_mean_pi, Y = Y,
                           n_lags = n_lags, nu = prior_nu)
  prior_pi <- priors$prior_pi
  prior_pi_omega <- priors$prior_pi_omega
  prior_s <- priors$prior_s

  # For the smoothing

  burn_in <-  gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                           Y, d, n_reps = n_burnin, n_fcst = NULL, Lambda, check_roots = TRUE,
                           d_fcst = NULL, init_Pi = t(prior_pi), init_psi = prior_mean_psi, smooth_state = FALSE)

  main_run <- gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_reps, n_fcst = 24, Lambda,
                            d_fcst = d_fcst,
                            init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                            init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],],
                            init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], ...)
  main_run$call <- fun_call
  return(main_run)

}

mfbvar2 <- function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_mean_pi,
                   lambda1, lambda2, prior_nu = NULL, prior_mean_psi, prior_var_psi, ...) {

  fun_call <- match.call()
  n_vars <- ncol(Y)
  n_T <- nrow(Y)
  n_T_ <- n_T - n_lags

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_pi_sigma(lambda1 = lambda1, lambda2 = lambda2, prior_mean = prior_mean_pi, Y = Y,
                           n_lags = n_lags, nu = prior_nu)
  prior_pi <- priors$prior_pi
  prior_pi_omega <- priors$prior_pi_omega
  prior_s <- priors$prior_s

  # For the smoothing
  M_Lambda <- build_M_Lambda(Y[-(1:n_lags), ], Lambda, n_vars, n_lags, n_T_)
  lH  <- M_Lambda

  burn_in <-  gibbs_sampler2(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_burnin, n_fcst = NULL, lH, check_roots = TRUE,
                            d_fcst = NULL, init_Pi = t(prior_pi), init_psi = prior_mean_psi, smooth_state = FALSE)

  main_run <- gibbs_sampler2(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_reps, n_fcst = 24, lH,
                            d_fcst = d_fcst,
                            init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                            init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],],
                            init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], ...)
  main_run$call <- fun_call
  return(main_run)

}

