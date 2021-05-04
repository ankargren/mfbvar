#' Get info about parameters
#'
#' Function that provides character vectors indicating different types of
#' parameters needed for a certain model
#'
#' @param minn logical for \code{minn} prior
#' @param ss logical for \code{ss} prior
#' @param ssng logical for \code{ssng} prior
#' @param dl logical for \code{dl} prior
#' @param iw logical for \code{iw} variance prior
#' @param diffuse logical for \code{diffuse} variance prior
#' @param csv logical for \code{csv} variance structure
#' @param fsv logical for \code{fsv} variance structure
#'
#' @return A list with components:
#' \item{prior_params}{Components in the prior object that need to be extracted}
#' \item{required_params}{Components that must be non-null}
#' \item{retrieved_params}{Components that are retrieved in initialization}
#' \item{params}{Parameters that are sampled}

get_params_info <- function(minn = FALSE,
                            ss = FALSE,
                            ssng = FALSE,
                            dl = FALSE,
                            iw = FALSE,
                            diffuse = FALSE,
                            csv = FALSE,
                            fsv = FALSE) {

  # Needed for all
  required_params <- c("Y", "n_lags", "n_burnin", "n_reps")

  retrieved_params <- c(
    "n_vars", "n_q", "T_b", "n_pseudolags",
    "n_T", "n_T_", "Z_1"
  )

  prior_params <- c(
    "freq", "verbose", "check_roots", "n_fcst", "n_thin",
    "freqs", "Lambda_", "lambda1", "lambda3",
    "lambda4", "prior_Pi_AR1", "block_exo"
  )

  params <- c("Z", "Pi")

  if (!fsv) {
    params <- c(params, "Sigma")
  }
  if (fsv) {
    required_params <- c(required_params, "n_fac")
    params <- c(
      params, "mu", "sigma", "phi", "facload", "f", "latent",
      "latent0"
    )
    prior_params <- c(
      prior_params, "priormu", "priorphiidi", "priorphifac",
      "priorsigmaidi", "priorsigmafac", "priorfacload",
      "restrict"
    )
  }
  if (csv) {
    params <- c(params, "phi", "latent", "latent0", "sigma")
    prior_params <- c(prior_params, "prior_phi", "prior_sigma2")
  }

  if (minn || dl) {
    prior_params <- c(prior_params, "lambda4")
  }

  if (dl) {
    params <- c(params, "slice")
  }

  if (ss || ssng) {
    prior_params <- c(prior_params, "d_fcst")
    required_params <- c(required_params, "prior_psi_mean", "d")
    params <- c(params, "psi")
    retrieved_params <- c(retrieved_params, "n_determ")
    if (!ssng) {
      required_params <- c(required_params, "prior_psi_Omega")
    }
  }
  if (ssng) {
    prior_params <- c(prior_params, "prior_ng", "s")
    params <- c(params, "omega", "lambda_mu", "phi_mu")
  }
  if (diffuse || fsv) {
    prior_params <- c(prior_params, "lambda2", "block_exo")
  }

  prior_params <- sort(c(required_params, prior_params))
  required_params <- sort(required_params)
  retrieved_params <- sort(retrieved_params)
  params <- sort(params)


  return(list(
    prior_params = prior_params,
    required_params = required_params,
    retrieved_params = retrieved_params,
    params = params
  ))
}
