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
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar n_fcst TRUE
#' @param ... Other arguments to pass to \code{\link{gibbs_sampler}}.
#' @template man_template
#' @details \code{mfbvar} calls \code{\link{gibbs_sampler}} (implemented in C++)

mfbvar <- function(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_Pi_AR1,
                   lambda1, lambda2, prior_nu = NULL, prior_psi_mean, prior_psi_Omega, n_fcst) {

  stopifnot(is.matrix(Y), is.matrix(d), is.matrix(prior_psi_Omega))
  stopifnot(is.vector(n_lags), is.vector(n_burnin), is.vector(n_reps), is.vector(lambda1), is.vector(lambda2), is.vector(n_fcst))
  stopifnot(nrow(Y) == nrow(d), ncol(Y) == length(prior_Pi_AR1), ncol(Y) == length(prior_psi_mean),
            ncol(Y) == sqrt(prod(dim(prior_psi_Omega))))

  if (!is.numeric(Lambda)) {
    if (is.character(Lambda) && length(Lambda) == ncol(Y)) {
      Lambda <- build_Lambda(Lambda, n_lags)
    } else {
      stop("Lambda must be either a matrix of correct size or a character vector with the aggregation use.")
    }
  }

  if (!is.null(n_fcst)) {
    if (nrow(d_fcst) != n_fcst) {
      stop(paste0("d_fcst has ", nrow(d_fcst), " rows, but n_fcst is ", n_fcst, "."))
    }
  }

  fun_call <- match.call()
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_Pi_Sigma(lambda1 = lambda1, lambda2 = lambda2, prior_Pi_AR1 = prior_Pi_AR1, Y = Y,
                           n_lags = n_lags, nu = prior_nu)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_s <- priors$prior_s

  # For the smoothing
  cat(paste0("############################################\nRunning the burn-in sampler with ", n_burnin, " draws\n"))
  start_burnin <- Sys.time()
  burn_in <-  gibbs_sampler(prior_Pi_mean, prior_Pi_Omega, prior_nu, prior_s, prior_psi_mean, prior_psi_Omega,
                           Y, d, n_reps = n_burnin, n_fcst = NULL, Lambda, check_roots = TRUE,
                           d_fcst = NULL, init_Pi = t(prior_Pi_mean), init_psi = prior_psi_mean, smooth_state = FALSE)
  end_burnin <- Sys.time()
  time_diff <- end_burnin - start_burnin
  cat(paste0("\nTime elapsed for drawing ", n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
             attr(time_diff, "units"), "\n"))
  cat(paste0("\nMoving on to ",
             n_reps, " replications in the main chain\n", ifelse(!is.null(n_fcst), paste0("Making forecasts ", n_fcst, " steps ahead"), NULL), "\n"))

  main_run <- gibbs_sampler(prior_Pi_mean, prior_Pi_Omega, prior_nu, prior_s, prior_psi_mean, prior_psi_Omega,
                            Y, d, n_reps = n_reps, n_fcst = n_fcst, Lambda,
                            d_fcst = d_fcst,
                            init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                            init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],],
                            init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]])
  main_run$call <- fun_call
  time_diff <- Sys.time() - start_burnin
  cat(paste0("\nTotal time elapsed: ", signif(time_diff, digits = 1), " ",
             attr(time_diff, "units")))
  return(main_run)

}

