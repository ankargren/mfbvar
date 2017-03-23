#' MF-BVAR with a steady-state prior
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @templateVar Y TRUE
#' @templateVar Lambda TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar lambda3 TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @templateVar verbose TRUE
#' @template man_template
#' @return An \code{mfbvar} object
#' @details \code{mfbvar} calls \code{\link{gibbs_sampler}} (implemented in C++)

mfbvar_schorf <- function(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_burnin, n_reps, verbose) {

  stopifnot(is.null(n_fcst) || is.vector(n_fcst))
  stopifnot(is.data.frame(Y) || is.matrix(Y), all(apply(Y, 2, is.numeric)))
  stopifnot(is.vector(n_lags), is.vector(n_burnin), is.vector(n_reps), is.vector(lambda1), is.vector(lambda2))
  stopifnot(ncol(Y) == length(prior_Pi_AR1))

  if (!is.numeric(Lambda)) {
    if (is.character(Lambda) && length(Lambda) == ncol(Y)) {
      Lambda <- build_Lambda(Lambda, n_lags)
    } else {
      stop("Lambda must be either a matrix of correct size or a character vector with the aggregation use.")
    }
  }

  if (!is.null(n_fcst)) {
    names_fcst <- paste0("fcst_", 1:n_fcst)
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

  original_Y <- Y
  Y <- as.matrix(Y)
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  # For the smoothing
  if (verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }
  burn_in <-  gibbs_sampler_schorf(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst = NULL, n_burnin,
                                 init_Pi = NULL, init_Sigma = NULL, init_Z = NULL,
                                 smooth_state = FALSE, check_roots = TRUE, verbose = TRUE)
  if (verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\n   Moving on to ",
               n_reps, " replications in the main chain\n", ifelse(!is.null(n_fcst), paste0("   Making forecasts ", n_fcst, " steps ahead"), " "), "\n\n"))
  }

  main_run <- gibbs_sampler_schorf(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_reps,
                                   init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                                   init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = FALSE, check_roots = TRUE, verbose)
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
  main_run$names_determ <- "const"
  main_run$n_burnin <- n_burnin
  main_run$prior_Pi_AR1 <- prior_Pi_AR1

  dimnames(main_run$Z) <- list(time = names_row,
                               variable = names_col,
                               iteration = 1:n_reps)
  dimnames(main_run$Pi) <- list(dep = names_col,
                                indep = c(paste0(rep(names_col, n_lags), ".l", rep(1:n_lags, each = n_vars)), "const"),
                                iteration = 1:n_reps)
  dimnames(main_run$Sigma) <- list(names_col,
                                   names_col,
                                   iteration = 1:n_reps)
    dimnames(main_run$psi) <- list(iteration = 1:n_reps,
                                 param = paste0(rep(names_col, 1), ".", rep(1, each = n_vars)))
  class(main_run) <- "mfbvar"
  return(main_run)

}


