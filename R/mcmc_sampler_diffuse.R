mcmc_sampler.mfbvar_minn_diffuse <- function(x, ...){

  required_params <- c("Y", "n_lags", "n_burnin", "n_reps")
  prior_params <- c(required_params,
                    "freq", "verbose", "check_roots", "n_fcst", "n_thin",
                    "freqs", "Lambda_", "lambda1", "lambda2", "lambda3",
                    "lambda4", "prior_Pi_AR1", "block_exo")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "Pi", "Sigma")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_dl_diffuse <- function(x, ...){

  required_params <- c("Y", "n_lags", "n_burnin", "n_reps")
  prior_params <- c(required_params,
                    "freq", "verbose", "check_roots", "n_fcst", "n_thin",
                    "freqs", "Lambda_", "lambda1", "lambda2", "lambda3",
                    "lambda4", "prior_Pi_AR1", "a", "gig", "n_cores",
                    "block_exo")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "Pi", "Sigma", "global", "aux", "local", "slice")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, dl = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ss_diffuse <- function(x, ...){

  required_params <- c("Y", "d", "n_lags", "n_burnin", "n_reps", "prior_psi_mean",
                       "prior_psi_Omega")
  prior_params <- c("Y", "d_fcst", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                    "lambda1", "lambda2", "lambda3", "prior_Pi_AR1",
                    "block_exo")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "psi", "Pi", "Sigma")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ss = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ssng_diffuse <- function(x, ...){

  required_params <- c("Y", "d", "n_lags", "n_burnin", "n_reps", "prior_psi_mean",
                       "prior_psi_Omega")
  prior_params <- c("Y", "d_fcst", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                    "lambda1", "lambda2", "lambda3", "prior_Pi_AR1",
                    "block_exo", "prior_ng", "s")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "psi", "Pi", "Sigma", "omega", "phi_mu", "lambda_mu")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ss = TRUE, ssng = TRUE, ...)

  return(out)
}
