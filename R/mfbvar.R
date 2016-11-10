mfbvar < function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, n_chains, prior_mean_pi, lambda1, lambda2, nu = NULL, prior_mean_psi, prior_var_psi, check_roots = TRUE) {

  n_vars <- ncol(Y)

  # Set priors
  if (is.null(nu)) {
    nu <- n_vars + 2
  }
  priors <- prior_pi_sigma(lambda1 = 0.2, lambda2 = 1, prior_mean = prior_mean_pi, Y = y,
                           n_lags = n_lags, nu = nu)
  prior_pi <- priors$prior_pi
  prior_pi_omega <- priors$prior_pi_omega
  prior_s <- priors$prior_s

  # For the smoothing
  M_Lambda <- build_M_Lambda(y[-(1:n_lags), ], Lambda, n_lags)
  lH  <- M_Lambda

  burn_in <- gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                        Y, d, n_reps = n_burnin, n_fcst = NULL, lH, check_roots = check_roots,
                        d_fcst = NULL)

  chain_init <- list(Pi_init = burn_in$Pi[,, n_burnin], Sigma_init = burn_in$Sigma[,, n_burnin],
                     psi_init = burn_in$psi[n_burnin, ], Z_init = burn_in$Z[,, n_burnin])

  cl <- makePSOCKcluster(n_chains)
  clusterSetRNGStream(cl, 2983468)
  registerDoParallel(cl)

  # Load the package
  clusterEvalQ(cl, {
    library(MFBVAR)
  })

  # Export the initializations from the burn in
  clusterExport(cl, chain_init, prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                Y, d, n_burnin, n_reps, n_fcst, lH, check_roots, d_fcst)

  main_run <- parLapply(cl, 1:n_chains, function(temp) {gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                                            Y, d, n_reps = round(n_burnin/5) + n_reps, n_fcst = n_fcst, lH, check_roots = check_roots,
                                            d_fcst = d_fcst)})




}
