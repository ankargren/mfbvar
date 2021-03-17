mcmc_sampler.mfbvar_minn_csv <- function(x, ...){

  params_info <- get_params_info(minn = TRUE, csv = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, csv = TRUE, ...)

  return(out)

}

mcmc_sampler.mfbvar_ss_csv <- function(x, ...) {

  params_info <- get_params_info(ss = TRUE, csv = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, ss = TRUE, csv = TRUE, ...)

  return(out)
}


mcmc_sampler.mfbvar_ssng_csv <- function(x, ...) {

  params_info <- get_params_info(ss = TRUE, ssng = TRUE, csv = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, ssng = TRUE, ss = TRUE, csv = TRUE, ...)

  return(out)
}


