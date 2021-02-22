mcmc_sampler.mfbvar_ss_iw <- function(x, ...) {

  required_params <- c("Y", "d", "prior_psi_mean", "prior_psi_Omega",
                       "n_lags", "n_burnin", "n_reps")
  prior_params <- c("Y", "freq", "verbose", "d", "d_fcst", "check_roots",
                    "n_lags", "n_fcst", "n_reps", "n_burnin", "n_thin",
                    "freqs", "Lambda_", "lambda1", "lambda3", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "n_determ", "Z_1")
  params <- c("Z", "psi", "Pi", "Sigma")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, ss = TRUE, iw = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ssng_iw <- function(x, ...) {

  required_params <- c("Y", "d", "prior_psi_mean",
                       "n_lags", "n_burnin", "n_reps")
  prior_params <- c("Y", "freq", "verbose", "d", "d_fcst", "check_roots",
                    "n_lags", "n_fcst", "n_reps", "n_burnin", "n_thin",
                    "freqs", "Lambda_", "lambda1", "lambda3",
                    "prior_Pi_AR1", "prior_ng", "s")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "n_determ", "Z_1")
  params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "Sigma")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, ssng = TRUE, ss = TRUE, iw = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_minn_iw <- function(x, ...){

  required_params <- c("Y", "n_lags", "n_burnin", "n_reps")
  prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "freqs", "Lambda_",
                    "lambda1", "lambda3", "lambda4", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "Pi", "Sigma")

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, iw = TRUE, ...)

  return(out)

}


