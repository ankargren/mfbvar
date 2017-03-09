#' MF-BVAR with a steady-state prior
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar d_fcst TRUE
#' @templateVar Lambda TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @templateVar prior_mean_Pi TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_mean_psi TRUE
#' @templateVar prior_var_psi TRUE
#' @param ... Other arguments to pass to \code{\link{gibbs_sampler}}.
#' @template man_template
#' @details \code{mfbvar} calls \code{\link{gibbs_sampler}} (implemented in C++)

mfbvar <- function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_mean_Pi,
                   lambda1, lambda2, prior_nu = NULL, prior_mean_psi, prior_var_psi, ...) {

  fun_call <- match.call()
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_Pi_Sigma(lambda1 = lambda1, lambda2 = lambda2, prior_mean = prior_mean_Pi, Y = Y,
                           n_lags = n_lags, nu = prior_nu)
  prior_Pi <- priors$prior_Pi
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_s <- priors$prior_s

  # For the smoothing

  burn_in <-  gibbs_sampler(prior_Pi, prior_Pi_Omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                           Y, d, n_reps = n_burnin, n_fcst = NULL, Lambda, check_roots = TRUE,
                           d_fcst = NULL, init_Pi = t(prior_Pi), init_psi = prior_mean_psi, smooth_state = FALSE)

  main_run <- gibbs_sampler(prior_Pi, prior_Pi_Omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_reps, n_fcst = 24, Lambda,
                            d_fcst = d_fcst,
                            init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                            init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],],
                            init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], ...)
  main_run$call <- fun_call
  return(main_run)

}


#' @template man_template
#' @inherit mfbvar
#' @details \code{mfbvar2} calls \code{\link{gibbs_sampler2}}
#'
mfbvar2 <- function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_mean_Pi,
                   lambda1, lambda2, prior_nu = NULL, prior_mean_psi, prior_var_psi, ...) {

  fun_call <- match.call()
  n_vars <- ncol(Y)
  n_T <- nrow(Y)
  n_T_ <- n_T - n_lags

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_Pi_Sigma(lambda1 = lambda1, lambda2 = lambda2, prior_mean = prior_mean_Pi, Y = Y,
                           n_lags = n_lags, nu = prior_nu)
  prior_Pi <- priors$prior_Pi
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_s <- priors$prior_s

  # For the smoothing
  M_Lambda <- build_M_Lambda(Y[-(1:n_lags), ], Lambda, n_vars, n_lags, n_T_)
  lH  <- M_Lambda

  burn_in <-  gibbs_sampler2(prior_Pi, prior_Pi_Omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_burnin, n_fcst = NULL, lH, check_roots = TRUE,
                            d_fcst = NULL, init_Pi = t(prior_Pi), init_psi = prior_mean_psi, smooth_state = FALSE)

  main_run <- gibbs_sampler2(prior_Pi, prior_Pi_Omega, prior_nu, prior_s, prior_mean_psi, prior_var_psi,
                            Y, d, n_reps = n_reps, n_fcst = 24, lH,
                            d_fcst = d_fcst,
                            init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                            init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],],
                            init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], ...)
  main_run$call <- fun_call
  return(main_run)

}

