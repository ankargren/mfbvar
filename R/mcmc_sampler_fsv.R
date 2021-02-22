mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){
  required_params <- c("Y", "n_lags", "n_burnin", "n_reps", "n_fac")
  prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                    "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                    "priorsigmafac", "priorfacload", "restrict",
                    "lambda1", "lambda2", "lambda3", "lambda4", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "Pi", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")
  out <- mfbvar_sampler(x, required_params, prior_params,
                               retrieved_params, params, fsv = TRUE, ...)

  return(out)
}


mcmc_sampler.mfbvar_dl_fsv <- function(x, ...){
  required_params <- c("Y", "n_lags", "n_burnin", "n_reps", "n_fac")
  prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                    "n_reps", "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                    "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                    "priorsigmafac", "priorfacload", "restrict",
                    "lambda1", "lambda2", "lambda3", "lambda4", "prior_Pi_AR1", "a",
                    "gig", "n_cores")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "Z_1")
  params <- c("Z", "Pi", "mu", "sigma", "phi", "facload", "f", "latent",
              "latent0", "global", "aux", "local", "slice")
  out <- mfbvar_sampler(x, required_params, prior_params,
                        retrieved_params, params, dl = TRUE, fsv = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ss_fsv <- function(x, ...){
  required_params <- c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags",
                       "n_burnin", "n_reps", "n_fac")
  prior_params <- c("Y", "freq", "verbose", "prior_psi_mean", "prior_psi_Omega",
                    "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                    "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                    "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                    "priorsigmafac", "priorfacload", "restrict",
                    "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "n_determ", "Z_1")
  params <- c("Z", "psi", "Pi", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")

  out <- mfbvar_sampler(x, required_params, prior_params,
                  retrieved_params, params, ss = TRUE, fsv = TRUE, ...)
  return(out)
}

mcmc_sampler.mfbvar_ssng_fsv <- function(x, ...){

  required_params <- c("Y", "d", "prior_psi_mean", "n_lags",
                       "n_burnin", "n_reps", "n_fac")
  prior_params <- c("Y", "freq", "verbose", "prior_psi_mean",
                    "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                    "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                    "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                    "priorsigmafac", "priorfacload", "restrict", "prior_ng", "s",
                    "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
  retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                        "n_T", "n_T_", "n_determ", "Z_1")
  params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")

  out <- mfbvar_sampler(x, required_params, prior_params,
                        retrieved_params, params, ssng = TRUE, ss = TRUE,
                        fsv = TRUE, ...)

  return(out)
}
