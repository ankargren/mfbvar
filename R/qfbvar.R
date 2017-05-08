#' MF-BVAR with a steady-state prior
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar d_fcst TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @templateVar verbose TRUE
#' @template man_template
#' @return An \code{mfbvar} object
#' @details \code{mfbvar} calls \code{\link{gibbs_sampler}} (implemented in C++)

qfbvar <- function(Y, d, d_fcst, prior_Pi_AR1, lambda1, lambda2, prior_nu = NULL, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, verbose, ...) {

  stopifnot(is.null(n_fcst) || is.vector(n_fcst))
  stopifnot(is.matrix(d), is.matrix(prior_psi_Omega))
  stopifnot(is.data.frame(Y) || is.matrix(Y), all(apply(Y, 2, is.numeric)))
  stopifnot(is.vector(n_lags), is.vector(n_burnin), is.vector(n_reps), is.vector(lambda1), is.vector(lambda2))
  stopifnot(nrow(Y) == nrow(d), ncol(Y) == length(prior_Pi_AR1), ncol(Y) == length(prior_psi_mean),
            ncol(Y) == sqrt(prod(dim(prior_psi_Omega))))

  if (!is.null(n_fcst)) {
    if (nrow(d_fcst) != n_fcst) {
      stop(paste0("d_fcst has ", nrow(d_fcst), " rows, but n_fcst is ", n_fcst, "."))
    } else {
      if (is.null(rownames(d_fcst))) {
        names_fcst <- paste0("fcst_", 1:n_fcst)
      } else {
        names_fcst <- rownames(d_fcst)
      }

    }
  } else {
    names_fcst <- NULL
  }

  fun_call <- match.call()
  if (is.null(rownames(Y))) {
    names_row <- 1:nrow(Y)
  } else {
    names_row <- rownames(Y)
  }

  if (is.null(colnames(Y))) {
    names_col <- 1:col(Y)
  } else {
    names_col <- colnames(Y)
  }

  if (is.null(colnames(d))) {
    names_determ <- paste0("d", 1:ncol(d))
  } else {
    names_determ <- colnames(d)
  }

  original_Y <- Y
  Y <- as.matrix(Y)
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_Pi_Sigma(lambda1 = lambda1, lambda2 = lambda2, prior_Pi_AR1 = prior_Pi_AR1, Y = Y,
                           n_lags = n_lags, prior_nu = prior_nu)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  # For the smoothing
  if (verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }
  burn_in <-  gibbs_sampler_qf(Y = Y, d = d, d_fcst = NULL, prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = prior_Pi_Omega,
                            prior_S = prior_S, prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega,
                            n_fcst = NULL, n_reps = n_burnin, check_roots = TRUE, verbose)
  if (verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\n   Moving on to ",
               n_reps, " replications in the main chain\n", ifelse(!is.null(n_fcst), paste0("   Making forecasts ", n_fcst, " steps ahead"), " "), "\n\n"))
  }

  main_run <- gibbs_sampler_qf(Y = Y, d = d, d_fcst = d_fcst, prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = prior_Pi_Omega,
                            prior_S = prior_S, prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega,
                            n_fcst = n_fcst, n_reps = n_reps, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                            init_psi = burn_in$psi[dim(burn_in$psi)[1],], check_roots = TRUE, verbose)
  main_run$call <- fun_call
  if (verbose) {
    time_diff <- Sys.time() - start_burnin
    cat(paste0("\n   Total time elapsed: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
  }

  if (!is.null(n_fcst)) {
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(main_run$n_T-main_run$n_lags+1):main_run$n_T]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col

  }


  main_run$names_row <- names_row
  main_run$names_col <- names_col
  main_run$names_fcst <- names_fcst
  main_run$names_determ <- names_determ
  main_run$n_burnin <- n_burnin
  main_run$prior_Pi_AR1 <- prior_Pi_AR1

  dimnames(main_run$Z) <- list(time = names_row,
                               variable = names_col)
  dimnames(main_run$Pi) <- list(dep = names_col,
                                indep = paste0(rep(names_col, n_lags), ".l", rep(1:n_lags, each = n_vars)),
                                iteration = 1:n_reps)
  dimnames(main_run$Sigma) <- list(names_col,
                                   names_col,
                                   iteration = 1:n_reps)
  n_determ <- dim(d)[2]
  dimnames(main_run$psi) <- list(iteration = 1:n_reps,
                                 param = paste0(rep(names_col, n_determ), ".", rep(names_determ, each = n_vars)))
  class(main_run) <- "qfbvar"
  return(main_run)

}

