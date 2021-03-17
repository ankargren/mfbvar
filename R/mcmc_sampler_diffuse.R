mcmc_sampler.mfbvar_minn_diffuse <- function(x, ...){

  params_info <- get_params_info(minn = TRUE, diffuse = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_dl_diffuse <- function(x, ...){

  params_info <- get_params_info(dl = TRUE, diffuse = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, dl = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ss_diffuse <- function(x, ...){

  params_info <- get_params_info(ss = TRUE, diffuse = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ss = TRUE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ssng_diffuse <- function(x, ...){

  params_info <- get_params_info(ss = TRUE, ssng = TRUE, diffuse = TRUE)
  required_params <- params_info$required_params
  prior_params <- params_info$prior_params
  retrieved_params <- params_info$retrieved_params
  params <- params_info$params

  out <- mfbvar_sampler(x, required_params, prior_params, retrieved_params,
                        params, diffuse = TRUE, ss = TRUE, ssng = TRUE, ...)

  return(out)
}
