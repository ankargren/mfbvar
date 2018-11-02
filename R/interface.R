#' Set priors for an mfbvar model
#'
#' Create an object storing all information needed for estimation, including data as well as model and prior specifications for both a Minnesota or steady-state prior.
#'
#' @templateVar Y TRUE
#' @templateVar freq TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @param d (Steady state only) Either a matrix with same number of rows as \code{Y} and \code{n_determ} number of columns containing the deterministic terms or a string \code{"intercept"} for requesting an intercept as the only deterministic
#' term.
#' @param d_fcst (Steady state only) The deterministic terms for the forecasting period (not used if \code{d = "intercept"}).
#' @param prior_psi_mean (Steady state only) Vector of length \code{n_determ*n_vars} with the prior means of the steady-state parameters.
#' @param prior_psi_Omega (Steady state only) Matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)} with the prior covariance of the steady-state parameters.
#' @param lambda3 (Minnesota only) Prior variance of the intercept.
#' @templateVar verbose TRUE
#' @templateVar smooth_state TRUE
#' @templateVar check_roots TRUE
#' @template man_template
#' @details The first arguments (\code{Y} through \code{n_reps}) must be set for the model to be estimated irrespective of the choice
#' of prior, but some have default values (which will produce warnings if relied upon).
#'
#' For the Minnesota prior, \code{lambda3} must also be set, but it too has a default that it relies on if not specified.
#'
#' For the steady-state prior, the deterministic matrix needs to be supplied, or a string indicating that the intercept should be
#' the only deterministic term. If the latter, also \code{d_fcst} is set to be intercept only. Otherwise, if forecasts are requested
#' (\code{n_fcst > 0}) also \code{d_fcst} needs to be provided. Finally, the prior moments for the steady-state parameters must also be
#' provided.
#'
#' The steady-state prior involves inverting the lag polynomial. For this reason, draws in which the largest eigenvalue
#' (in absolute value) of the lag polynomial is greater than 1 are discarded and new draws are made. The maximum number of
#' attempts is 1,000. The components in the output named \code{roots} and \code{num_tries} contain the largest roots and the
#' number of attempts, respectively, if \code{check_roots = TRUE} (the default).
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 100, n_reps = 100)
#' prior_obj <- update_prior(prior_obj, n_fcst = 4)
#' @seealso \code{\link{interval_to_moments}}, \code{\link{print.mfbvar_prior}}, \code{\link{summary.mfbvar_prior}}, \code{\link{estimate_mfbvar}}
set_prior <- function(Y, freq, prior_Pi_AR1 = rep(0, ncol(Y)), lambda1 = 0.2,
                      lambda2 = 1, n_lags, n_fcst = 0, n_burnin, n_reps,
                      d = NULL, d_fcst = NULL, prior_psi_mean = NULL,
                      prior_psi_Omega = NULL, lambda3 = 10000, verbose = FALSE,
                      smooth_state = FALSE, check_roots = TRUE) {
  prior_call <- mget(names(formals()),sys.frame(sys.nframe()))
  prior_call$supplied_args <- names(as.list(match.call()))[-1]
  ret <- check_prior(prior_call)
  class(ret) <- "mfbvar_prior"
  return(ret)
}
#' @rdname set_prior
#'
#' @param prior_obj an object of class \code{mfbvar_prior}
#' @param ... named arguments for prior attributes to update
update_prior <- function(prior_obj, ...) {
  if(!inherits(prior_obj, "mfbvar_prior")) {
    stop("The object to be updated is not of class mfbvar_prior.")
  }

  tmp <- list(...)
  nam <- names(tmp)
  for(iter in 1:length(tmp)) {
    eval(parse(text=paste0("prior_obj$",nam[iter]," = tmp[[iter]]")))
  }
  prior_obj$supplied_args <- union(prior_obj$supplied_args, names(as.list(match.call()))[-1])
  prior_obj <- check_prior(prior_obj)

  return(prior_obj)
}


#' @rdname set_prior
check_prior <- function(prior_obj) {
  if (!is.matrix(prior_obj$Y)) {
    if (!is.data.frame(prior_obj$Y)) {
      stop(paste0("Y is of class ", class(prior_obj$Y), ", but must be matrix or data frame."))
    } else {
      prior_obj$Y <- as.matrix(prior_obj$Y)
    }
  }

  if (nrow(prior_obj$Y) > max(unlist(apply(prior_obj$Y, 2, function(x) Position(is.na, x, nomatch = nrow(prior_obj$Y)))))) {
    stop("Y: remove final rows containing only NAs.")
  }

  intercept_flag <- FALSE

  if ("d" %in% prior_obj$supplied_args) {
    if (is.atomic(prior_obj$d) && length(prior_obj$d) == 1) {
      intercept_flag <- TRUE
      prior_obj$d <- matrix(1, nrow = nrow(prior_obj$Y), 1)
    } else if (all(prior_obj$d == 1)) {
      intercept_flag <- TRUE
    }

    prior_obj$intercept_flag <- intercept_flag

    if (!intercept_flag) {
      if ("d_fcst" %in% prior_obj$supplied_args) {
        if (!is.matrix(prior_obj$d_fcst)) {
          stop("d_fcst is of class", class(prior_obj$d_fcst), ", but must be a matrix.")
        }

        if (ncol(prior_obj$d) != ncol(prior_obj$d_fcst)) {
          stop("d has", ncol(prior_obj$d), " columns and d_fcst", ncol(prior_obj$d_fcst), ".")
        }
      }
    }

    if (nrow(prior_obj$Y) != nrow(prior_obj$d)) {
      stop("Y has ", nrow(prior_obj$Y), "rows and d ", nrow(prior_obj$d), "rows, but they must be equal.")
    }

    if (!is.matrix(prior_obj$d)) {
      stop("d is of class", class(prior_obj$d), ", but must be a matrix.")
    }
  }


  if ("freq" %in% prior_obj$supplied_args) {
    if (!(is.atomic(prior_obj$freq) && is.character(prior_obj$freq))) {
      stop("freq is of class ", class(prior_obj$freq), ", but it must be a character vector.")
    } else if (!all(prior_obj$freq %in% c("m", "q"))) {
      stop("Elements of freq must be 'm' or 'q'.")
    } else if (length(prior_obj$freq) != ncol(prior_obj$Y)) {
      stop("The length of freq is ", length(prior_obj$freq), ", but Y has ", ncol(prior_obj$Y), " columns.")
    } else if (which.max(prior_obj$freq == "m") > which.max(prior_obj$freq == "q")) {
      stop("Monthly variables must be placed before quarterly variables.")
    }
  } else {
    stop("freq: must be supplied.")
  }

  if (min(unlist(apply(prior_obj$Y[, prior_obj$freq == "m", drop = FALSE], 2, function(x) Position(is.na, x, nomatch = 9999999999)))) == 1) {
    stop("Y: monthly variables are NA at the beginning of the sample.")
  }

  if ("prior_Pi_AR1" %in% prior_obj$supplied_args) {
    if (is.atomic(prior_obj$prior_Pi_AR1)) {
      if (length(prior_obj$prior_Pi_AR1) == 1) {
        warning("prior_Pi_AR1: Recycling ", prior_obj$prior_Pi_AR1, " to use as prior mean for AR(1) coefficients for all variables.\n", call. = FALSE)
        prior_obj$prior_Pi_AR1 <- rep(prior_obj$prior_Pi_AR1, ncol(prior_obj$Y))
      }
    } else {
      stop("prior_Pi_AR1 must be a vector, but is now of class ", class(prior_obj$prior_Pi_AR1), ".")
    }
  } else {
    warning("prior_Pi_AR1: 0 used as prior mean for AR(1) coefficients.\n", call. = FALSE)
    prior_obj$prior_Pi_AR1 <- rep(0, ncol(prior_obj$Y))
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "prior_Pi_AR1")
  }


  if ("lambda1" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda1) || length(prior_obj$lambda1) > 1) {
      stop("lambda1 must be a vector with a single element.")
    }
  } else {
    warning("lambda1: ", prior_obj$lambda1, " used as the value for the overall tightness hyperparameter.\n", call. = FALSE)
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda1")
  }

  if ("lambda2" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda2) || length(prior_obj$lambda2) > 1) {
      stop("lambda2 must be a vector with a single element.")
    }
  } else {
    warning("lambda2: ", prior_obj$lambda2, " used as the value for the lag decay hyperparameter.\n", call. = FALSE)
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda2")
  }

  if ("lambda3" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda3) || length(prior_obj$lambda3) > 1) {
      stop("lambda3 must be a vector with a single element.")
    }
  } else {
    warning("lambda3: ", prior_obj$lambda3, " used for the constant's prior variance.\n", call. = FALSE)
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda3")
  }


  if ("prior_psi_mean" %in% prior_obj$supplied_args) {
    if (!(is.atomic(prior_obj$prior_psi_mean) || is.matrix(prior_obj$prior_psi_mean))) {
      stop("prior_psi_mean must be a vector or matrix with one row or column.")
    }
    if (is.atomic(prior_obj$prior_psi_mean)) {
      if (length(prior_obj$prior_psi_mean) != ncol(prior_obj$Y)) {
        stop("prior_psi_mean has ", length(prior_obj$prior_psi_mean), " elements, but there are ", ncol(prior_obj$Y), " variables in Y.")
      }
    }
    if (is.matrix(prior_obj$prior_psi_mean)) {
      if (!any(dim(prior_obj$prior_psi_mean) == 1)) {
        stop("prior_psi_mean must be a vector or matrix with one row or column.")
      } else {
        prior_obj$prior_psi_mean <- c(prior_obj$prior_psi_mean)
      }
    }
  }


  if ("prior_psi_Omega" %in% prior_obj$supplied_args) {
    if (!is.matrix(prior_obj$prior_psi_Omega)) {
      stop("prior_psi_Omega must be a matrix.")
    } else {
      if (dim(prior_obj$prior_psi_Omega)[1] != dim(prior_obj$prior_psi_Omega)[2]) {
        stop("prior_psi_Omega must be a positive-definite symmetric matrix.")
      }
    }
  }


  if ("n_lags" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_lags) || length(prior_obj$n_lags) > 1) {
      stop("n_lags must be a vector with a single element.")
    }
  } else {
    stop("n_lags: No lag length specified.\n")
  }

  if ("n_fcst" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_fcst) || length(prior_obj$n_fcst) > 1) {
      stop("n_fcst must be a vector with a single element.")
    } else {
      if (intercept_flag) {
        prior_obj$d_fcst <- matrix(1, nrow = prior_obj$n_fcst, 1)
      }
    }
  } else {
    warning("n_fcst: ", prior_obj$n_fcst, " used for the number of forecasts to compute.\n", call. = FALSE)
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "n_fcst")
  }



  if ("n_burnin" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_burnin) || length(prior_obj$n_burnin) > 1) {
      stop("n_burnin must be a vector with a single element.")
    }
  } else {
    stop("n_burnin: Number of burn-in draws to use not specified.\n")
  }

  if ("n_reps" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_reps) || length(prior_obj$n_reps) > 1) {
      stop("n_reps must be a vector with a single element.")
    }
  } else {
    stop("n_reps: Number of draws to use in main chain not specified.\n")
  }

  if (!is.logical(prior_obj$smooth_state)) {
    stop("smooth_state: must be logical.\n")
  }

  if (!is.logical(prior_obj$check_roots)) {
    stop("check_roots: must be logical.\n")
  }

  return(prior_obj)
}

#' Print method for mfbvar_prior
#'
#' Printing method for object of class mfbvar_prior, checking if information
#' in the prior is sufficient for estimating models.
#' @param x prior object (class \code{mfbvar_prior})
#' @param ... additional arguments (currently unused)
#' @details The print method checks whether the steady-state and Minnesota
#'   priors can be run with the current specification. This check is minimal in
#'   the sense that it checks only prior elements with no defaults, and it only
#'   checks for estimation and not forecasting (for which the steady-state prior
#'   requires additional information).
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{estimate_mfbvar}}, \code{\link{summary.mfbvar_prior}}
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 100, n_reps = 100)
#' print(prior_obj)
print.mfbvar_prior <- function(x, ...) {
  cat("The following elements of the prior have not been set: \n", names(sapply(x, is.null))[sapply(x, is.null)])
  cat("\n\n")
  cat("Checking if steady-state prior can be run... ")
  if (!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$prior_psi_Omega) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n")
  }

  cat("Checking if Minnesota prior can be run... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
  }
}

#' Summary method for mfbvar_prior
#'
#' summary method for object of class mfbvar_prior, showing some basic
#' information regarding the contents of the prior.
#' @param object prior object (class \code{mfbvar_prior})
#' @param ... additional arguments (currently unused)
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{estimate_mfbvar}}, \code{\link{print.mfbvar_prior}}
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 100, n_reps = 100)
#' summary(prior_obj)
summary.mfbvar_prior <- function(object, ...) {
  cat("PRIOR SUMMARY\n")
  cat("----------------------------\n")
  cat("Required elements:\n")
  cat("  Y:", ncol(object$Y), "variables,", nrow(object$Y), "time points\n")
  cat("  freq:", object$freq, "\n")
  cat("  prior_Pi_AR1:", object$prior_Pi_AR1, "\n")
  cat("  lambda1:", object$lambda1, "\n")
  cat("  lambda2:", object$lambda2, "\n")
  cat("  n_lags:", object$n_lags, "\n")
  cat("  n_fcst:", object$n_fcst, "\n")
  cat("  n_burnin:", object$n_burnin, "\n")
  cat("  n_reps:", object$n_reps, "\n")
  cat("----------------------------\n")
  cat("Steady-state-specific elements:\n")
  cat("  d:", ifelse(is.null(object$d), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(ncol(object$d), "deterministic variables"))),"\n")
  cat("  d_fcst:", ifelse(is.null(object$d_fcst), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(nrow(object$d_fcst), "forecasts, ", ncol(object$d), "deterministic variables"))),"\n")
  cat("  prior_psi_mean:", ifelse(is.null(object$prior_psi_mean), "<missing>", "vector of prior steady-state means"), "\n")
  cat("  prior_psi_Omega:", ifelse(is.null(object$prior_psi_Omega), "<missing>", "prior covariance matrix of prior steady states"), "\n")
  cat("----------------------------\n")
  cat("Minnesota-specific elements:\n")
  cat("  lambda3:", ifelse(is.null(object$lambda3), "<missing> (will rely on default)", object$lambda3), "\n")
  cat("----------------------------\n")
  cat("Other:\n")
  cat("  verbose:", object$verbose, "\n")
  cat("  smooth_state:", object$smooth_state, "\n")
  cat("  check_roots:", object$check_roots, "\n")

}


#' Mixed-frequency Bayesian VAR
#'
#' The main function for estimating a mixed-frequency BVAR.
#'
#' @param mfbvar_prior a \code{mfbvar_prior} object
#' @param prior_type either \code{"ss"} (steady-state prior) or \code{"minn"} (Minnesota prior)
#' @param ... additional arguments to \code{update_prior} (if \code{mfbvar_prior} is \code{NULL}, the arguments are passed on to \code{set_prior})
#' @return an object of class \code{mfbvar} and \code{mfbvar_ss}/\code{mfbvar_minn} containing posterior quantities as well as the prior object
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{predict.mfbvar}}, \code{\link{plot.mfbvar_minn}},
#' \code{\link{plot.mfbvar_ss}}, \code{\link{summary.mfbvar_minn}}, \code{\link{summary.mfbvar_ss}}
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' @return The prior values used are carried forward and returned with \code{NULL} if not used/existing. New components are:
#' \item{Pi}{Array of dynamic coefficient matrices (\eqn{\Pi}) from the main chain; \code{Pi[,, r]} is the \code{r}th draw}
#' \item{Sigma}{Array of covariance matrices (\eqn{\Sigma}) from the main chain; \code{Sigma[,, r]} is the \code{r}th draw}
#' \item{psi}{Matrix of steady-state parameter vectors (\eqn{\psi}) from the main chain; \code{psi[r,]} is the \code{r}th draw}
#' \item{Z}{Array of monthly process (\eqn{Z}) matrices from the main chain; \code{Z[,, r]} is the \code{r}th draw}
#' \item{roots}{The maximum eigenvalue of the lag polynomial (if \code{check_roots = TRUE})}
#' \item{num_tries}{The number of attempts for drawing a stationary \eqn{\Pi} (only relevant if \code{prior_type = "ss"})}
#' \item{Z_fcst}{Array of monthly forecasts from the main chain; \code{Z_fcst[,, r]} is the \code{r}th forecast. The first \code{n_lags}
#' rows are taken from the data to offer a bridge between observations and forecasts and for computing nowcasts (i.e. with ragged edges).}
#' \item{smoothed_Z}{The smoothed estimates (if \code{smooth_state = TRUE})}
#' @references
#' Schorfheide, F., & Song, D. (2015) Real-Time Forecasting With a Mixed-Frequency VAR. \emph{Journal of Business & Economic Statistics}, 33(3), 366--380. \url{http://dx.doi.org/10.1080/07350015.2014.954707}\cr
#' Ankargren, S., Unosson, M., & Yang, Y. (2018) A Mixed-Frequency Bayesian Vector Autoregression with a Steady-State Prior. Working Paper, Department of Statistics, Uppsala University No. 2018:3.
estimate_mfbvar <- function(mfbvar_prior = NULL, prior_type, ...) {
  if (hasArg(mfbvar_prior)) {
    if (!inherits(mfbvar_prior, "mfbvar_prior")) {
      stop("mfbvar_prior must be of class mfbvar_prior.")
    } else {
      if (length(list(...)) > 0) {
        mfbvar_prior <- update_prior(mfbvar_prior, ...)
      }
    }
  } else {
    mfbvar_prior <- set_prior(...)
  }

  if (!(prior_type %in% c("ss", "minn"))) {
    stop("prior_type must be 'ss' or 'minn'.")
  }

  class(mfbvar_prior) <- c(class(mfbvar_prior), paste0("mfbvar_", prior_type))

  if (mfbvar_prior$verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", mfbvar_prior$n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }

  burn_in <- mcmc_sampler(update_prior(mfbvar_prior, n_fcst = 0, smooth_state = FALSE), n_reps = mfbvar_prior$n_burnin)

  if (mfbvar_prior$verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", mfbvar_prior$n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\n   Moving on to the main chain with ",
               mfbvar_prior$n_reps, " draws \n\n", ifelse(mfbvar_prior$n_fcst > 0, paste0("   Making forecasts ", mfbvar_prior$n_fcst, " steps ahead"), ""), "\n\n"))
  }

  main_run <- mcmc_sampler(mfbvar_prior, n_reps = mfbvar_prior$n_reps+1, init = burn_in$init)

  if (mfbvar_prior$verbose) {
    time_diff <- Sys.time() - start_burnin
    cat(paste0("\n   Total time elapsed: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
  }


  n_vars <- ncol(mfbvar_prior$Y)
  if (is.null(rownames(mfbvar_prior$Y))) {
    names_row <- 1:nrow(mfbvar_prior$Y)
  } else {
    names_row <- rownames(mfbvar_prior$Y)
  }

  if (is.null(colnames(mfbvar_prior$Y))) {
    names_col <- 1:ncol(mfbvar_prior$Y)
  } else {
    names_col <- colnames(mfbvar_prior$Y)
  }
  if (mfbvar_prior$n_fcst > 0) {
    names_fcst <- paste0("fcst_", 1:mfbvar_prior$n_fcst)
    main_run$Z_fcst <- main_run$Z_fcst[,,-1]
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(main_run$n_T-main_run$n_lags+1):main_run$n_T]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col
  } else {
    names_fcst <- NULL
  }


  main_run$names_row <- names_row
  main_run$names_col <- names_col
  main_run$names_fcst <- names_fcst
  main_run$mfbvar_prior <- mfbvar_prior
  main_run$Pi <- main_run$Pi[,,-1]
  main_run$Sigma <- main_run$Sigma[,,-1]
  main_run$Z <- main_run$Z[,,-1]

  dimnames(main_run$Z) <- list(time = names_row,
                               variable = names_col,
                               iteration = 1:mfbvar_prior$n_reps)
  dimnames(main_run$Sigma) <- list(names_col,
                                   names_col,
                                   iteration = 1:mfbvar_prior$n_reps)

  if (prior_type == "ss") {
    if (is.null(colnames(mfbvar_prior$d))) {
      names_determ <- paste0("d", 1:ncol(mfbvar_prior$d))
    } else {
      names_determ <- colnames(mfbvar_prior$d)
    }
    main_run$names_determ <- names_determ
    n_determ <- dim(mfbvar_prior$d)[2]
    main_run$psi <- main_run$psi[-1, ]
    dimnames(main_run$psi) <- list(iteration = 1:mfbvar_prior$n_reps,
                                   param = paste0(rep(names_col, n_determ), ".", rep(names_determ, each = n_vars)))
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars)),
                                  iteration = 1:mfbvar_prior$n_reps)
  } else {
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = c(paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars)), "const"),
                                  iteration = 1:mfbvar_prior$n_reps)
  }

  class(main_run) <- c("mfbvar", paste0("mfbvar_", prior_type))
  return(main_run)
}


#' Printing method for class mfbvar_ss
#'
#' Method for printing \code{mfbvar_ss} objects.
#'
#' @param x object of class \code{mfbvar_ss}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20)
#' prior_intervals <- matrix(c(-0.1, 0.1,
#'                             0.4, 0.6), ncol = 2, byrow = TRUE)
#' psi_moments <- interval_to_moments(prior_intervals)
#' prior_psi_mean <- psi_moments$prior_psi_mean
#' prior_psi_Omega <- psi_moments$prior_psi_Omega
#' prior_obj <- update_prior(prior_obj,
#'                           prior_psi_mean = prior_psi_mean,
#'                           prior_psi_Omega = prior_psi_Omega)
#' mod_ss <- estimate_mfbvar(prior_obj, prior_type = "ss")
#' print(mod_ss)

print.mfbvar_ss <- function(x, ...){
  cat(paste0("Mixed-frequency steady-state BVAR with:\n", ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "), "\n", x$n_lags, " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecasted\n",
             x$n_reps, " draws used in main chain"))
}

#' Printing method for class mfbvar_minn
#'
#' Method for printing \code{mfbvar_minn} objects.
#'
#' @param x object of class \code{mfbvar_minn}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' print(mod_minn)

print.mfbvar_minn <- function(x, ...){
  cat(paste0("Mixed-frequency Minnesota BVAR with:\n", ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "), "\n", x$n_lags, " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecasted\n",
             x$n_reps, " draws used in main chain"))
}

#' Plotting method for class \code{mfbvar_ss}
#'
#' Method for plotting \code{mfbvar_ss} objects.
#' @param x object of class \code{mfbvar_ss}
#' @param plot_start Time period (number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param ss_level A vector with the lower and upper quantiles for the posterior steady-state intervals.
#' @param pred_level A vector with the lower and upper quantiles for the forecast intervals.
#' @param nrow_facet an integer giving the number of rows to use in the facet
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20,
#'                        n_fcst = 4)
#'
#' prior_intervals <- matrix(c(-0.1, 0.1,
#'                             0.4, 0.6), ncol = 2, byrow = TRUE)
#' psi_moments <- interval_to_moments(prior_intervals)
#' prior_psi_mean <- psi_moments$prior_psi_mean
#' prior_psi_Omega <- psi_moments$prior_psi_Omega
#' prior_obj <- update_prior(prior_obj,
#'                           prior_psi_mean = prior_psi_mean,
#'                           prior_psi_Omega = prior_psi_Omega)
#'
#' mod_ss <- estimate_mfbvar(prior_obj, prior_type = "ss")
#' plot(mod_ss)

plot.mfbvar_ss <- function(x, plot_start = NULL, ss_level = c(0.025, 0.975),
                        pred_level = c(0.10, 0.90), nrow_facet = NULL, ...){
  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(x$n_T-x$n_fcst*5, 0):x$n_T
    }  else {
      plot_range <- 1:x$n_T
    }
  } else {
    plot_range <- plot_start:x$n_T
  }

  if (is.null(ss_level)) {
    ss_level <- c(0.025, 0.975)
  }
  if (is.null(pred_level)) {
    pred_level <- c(0.10, 0.90)
  }

  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  names_row <- if (is.null(x$names_row)) 1:x$n_T else x$names_row
  p <- ggplot(mapping = aes(x = time))

  if (x$n_fcst > 0) {
    ss_lower  <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
    ss_median <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
    ss_upper  <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))
  } else {

    ss_lower  <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
    ss_median <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
    ss_upper  <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))

  }

  if (!is.null(x$psi)) {
    ss <- data.frame(expand.grid(time = plot_range[1]:(plot_range[length(plot_range)]+x$n_fcst), names_col = names_col), lower = c(ss_lower), median = c(ss_median),
                     upper = c(ss_upper))
    ss$value <- c(rbind(as.matrix(x$Y[plot_range,]), matrix(NA, nrow = x$n_fcst, ncol = x$n_vars)))
    ss_excl <- c(rbind(!is.na(as.matrix(x$Y[plot_range,])), matrix(TRUE, nrow = x$n_fcst, ncol = x$n_vars)))
    #ss <- na.omit(ss)

    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "#bdbdbd"),  alpha =1) +
      geom_line(data = na.omit(ss), aes(y = value))
  }
  if (x$n_fcst > 0) {
    preds <- predict(x, pred_quantiles = c(pred_level[1], 0.5, pred_level[2]))
    fcst <- data.frame(expand.grid(time = (x$n_T+x$n_fcst-nrow(preds[[1]])+1):(x$n_T+x$n_fcst), names_col = x$names_col),
                       lower = c(preds[[1]]), median = c(preds[[2]]),
                       upper = c(preds[[3]]))
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(time = (1:nrow(x$Y))[last_pos[i]],
                                     names_col = names_col[i],
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]))
    }
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "grey90"), linetype = "dotted", color = "black",
                         alpha = 0.75) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = "Note: The forecasts are for the underlying variable.")
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }

  if (x$n_fcst > 0) {
    p <- p + scale_fill_manual(values = c("grey90" = "grey90", "#bdbdbd" = "#bdbdbd"),
                               label = c("grey90" = paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"),
                               "#bdbdbd" = paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Forecasts and posterior steady state intervals",
           y = "Value",
           x = "Time") +
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd", "grey90"),
                                                     linetype = c("blank", "dotted"))))
  } else {
    p <- p + scale_fill_manual(values = c("#bdbdbd"),
                               label = c(paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Steady state intervals",
           y = "Value",
           x = "Time")+
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd"),
                                                     linetype = c("blank"))))
  }


  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~names_col, scales = "free_y")
  } else {
    p <- p + facet_wrap(~names_col, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position="bottom")
  breaks <- ggplot_build(p)$layout$panel_ranges[[1]]$x.major_source
  if (length(which(!(breaks %in% 1:x$n_T))) > 0) {
    breaks <- breaks[-which(!(breaks %in% 1:x$n_T))]
  }
  p + scale_x_continuous(breaks = breaks,
                         labels = names_row[breaks])
}

#' Plotting method for class \code{mfbvar_minn}
#'
#' Method for plotting mfbvar objects.
#' @param x object of class \code{mfbvar_minn}
#' @param plot_start Time period (number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param pred_level A vector with the lower and upper quantiles for the forecast intervals.
#' @param nrow_facet an integer giving the number of rows to use in the facet
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20, n_fcst = 4)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' plot(mod_minn)

plot.mfbvar_minn <- function(x, plot_start = NULL, pred_level = c(0.10, 0.90), nrow_facet = NULL, ...){
  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(x$n_T-x$n_fcst*5, 0):x$n_T
    }  else {
      plot_range <- 1:x$n_T
    }
  } else {
    plot_range <- plot_start:x$n_T
  }

  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  names_row <- if (is.null(x$names_row)) 1:x$n_T else x$names_row
  p <- ggplot(mapping = aes(x = time))

  ss <- data.frame(expand.grid(time = plot_range, names_col = names_col))
  ss$value <- c(as.matrix(x$Y[plot_range,]))
  ss <- na.omit(ss)

  p <- p +
    geom_line(data = ss, aes(y = value))

  if (x$n_fcst > 0) {
    preds <- predict(x, pred_quantiles = c(pred_level[1], 0.5, pred_level[2]))
    fcst <- data.frame(expand.grid(time = (x$n_T+x$n_fcst-nrow(preds[[1]])+1):(x$n_T+x$n_fcst), names_col = x$names_col),
                       lower = c(preds[[1]]), median = c(preds[[2]]),
                       upper = c(preds[[3]]))
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(time = (1:nrow(x$Y))[last_pos[i]],
                                     names_col = names_col[i],
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]))
    }
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <-p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                           fill = "grey90"), linetype = "dotted", color = "black",
                          alpha = 0.75) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = "Note: The forecasts are for the underlying variable.")
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }



  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~names_col, scales = "free_y")
  } else {
    p <- p + facet_wrap(~names_col, scales = "free_y", nrow = nrow_facet)
  }
  p <- p  +
    scale_fill_manual(values = c("grey90" = "grey90"),
                      label = c(paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"))) +
    labs(fill = "Intervals",
         title = "Forecast intervals",
         y = "Value",
         x = "Time") +
    theme_minimal() +
    theme(legend.position="bottom")
  breaks <- ggplot_build(p)$layout$panel_ranges[[1]]$x.major_source
  if (length(which(!(breaks %in% 1:x$n_T))) > 0) {
    breaks <- breaks[-which(!(breaks %in% 1:x$n_T))]
  }
  p + scale_x_continuous(breaks = breaks,
                         labels = names_row[breaks])
}

#' Summary method for class \code{mfbvar_ss}
#'
#' Method for summarizing \code{mfbvar_ss} objects.
#'
#' @param object object of class \code{mfbvar_ss}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20,
#'                        n_fcst = 4)
#'
#' prior_intervals <- matrix(c(-0.1, 0.1,
#'                             0.4, 0.6), ncol = 2, byrow = TRUE)
#' psi_moments <- interval_to_moments(prior_intervals)
#' prior_psi_mean <- psi_moments$prior_psi_mean
#' prior_psi_Omega <- psi_moments$prior_psi_Omega
#' prior_obj <- update_prior(prior_obj,
#'                           prior_psi_mean = prior_psi_mean,
#'                           prior_psi_Omega = prior_psi_Omega)
#'
#' mod_ss <- estimate_mfbvar(prior_obj, prior_type = "ss")
#' summary(mod_ss)
#'
summary.mfbvar_ss <- function(object, ...) {
  post_Pi <- apply(object$Pi, c(1, 2), mean)
  rownames(post_Pi) <- object$names_col
  colnames(post_Pi) <- paste0(rep(object$names_col, object$n_lags), ".", rep(1:object$n_lags, each = object$n_vars))
  post_Sigma <- apply(object$Sigma, c(1, 2), mean)
  rownames(post_Sigma) <- object$names_col
  colnames(post_Sigma) <- object$names_col
  post_psi <- matrix(colMeans(object$psi), ncol = object$n_determ)
  rownames(post_psi) <- object$names_col
  colnames(post_psi) <- object$names_determ
  print(object, ...)
  cat("\n\n#########################\nPosterior means computed\n\nPi:\n")
  print(post_Pi)
  cat("\n\nSigma:\n")
  print(post_Sigma)
  cat("\n\nPsi:\n")
  print(post_psi)
  ret_list <- list(post_Pi = post_Pi, post_Sigma = post_Sigma, post_Psi = post_psi)
}

#' Summary method for class \code{mfbvar_minn}
#'
#' Method for summarizing \code{mfbvar_minn} objects.
#'
#' @param object object of class \code{mfbvar_minn}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20, n_fcst = 4)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' summary(mod_minn)
#'
summary.mfbvar_minn <- function(object, ...) {
  post_Pi <- apply(object$Pi[,-ncol(object$Pi[,,1]),], c(1, 2), mean)
  rownames(post_Pi) <- object$names_col
  colnames(post_Pi) <- paste0(rep(object$names_col, object$n_lags), ".", rep(1:object$n_lags, each = object$n_vars))
  post_Sigma <- apply(object$Sigma, c(1, 2), mean)
  rownames(post_Sigma) <- object$names_col
  colnames(post_Sigma) <- object$names_col
  post_psi <- matrix(rowMeans(object$Pi[,ncol(object$Pi[,,1]),]), ncol = 1)
  rownames(post_psi) <- object$names_col
  colnames(post_psi) <- "const"
  print(object, ...)
  cat("\n\n#########################\nPosterior means computed\n\nPi:\n")
  print(post_Pi)
  cat("\n\nSigma:\n")
  print(post_Sigma)
  cat("\n\nIntercept:\n")
  print(post_psi)
  ret_list <- list(post_Pi = post_Pi, post_Sigma = post_Sigma, post_Psi = post_psi)
}

#' Predict method for class \code{mfbvar}
#'
#' Method for predicting \code{mfbvar} objects.
#'
#' @param object object of class mfbvar
#' @param pred_quantiles The quantiles of the posterior predictive distribution to use.
#' @param tidy If results should be tidy or not.
#' @param ... Currently not in use.
#' @template man_template
#' @details Note that this requires that forecasts were made in the original \code{mfbvar} call.
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20, n_fcst = 4)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' predict(mod_minn)
#' predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE)
predict.mfbvar <- function(object, pred_quantiles = c(0.10, 0.50, 0.90), tidy = FALSE, ...) {
  if (object$n_fcst==0) {
    stop("No forecasts exist in the provided object.")
  }

  final_non_na <- min(unlist(apply(object$Y, 2, function(x) Position(is.na, x, nomatch = nrow(object$Y))))[object$mfbvar_prior$freq == "m"])
  final_fcst <- object$n_lags- (nrow(object$Y)-final_non_na)
  if (final_fcst >= 1) {
    incl_fcst <- final_fcst:(object$n_lags + object$n_fcst)
  } else {
    incl_fcst <- 1:(object$n_lags + object$n_fcst)
  }
  if (tidy == FALSE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[incl_fcst,,], c(1, 2), quantile, prob = xx))
    names(ret_list) <- paste0("quantile_", pred_quantiles*100)
    return(ret_list)
  } else if (tidy == TRUE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[incl_fcst,,], c(1, 2), quantile, prob = xx))
    ret_tidy <- cbind(value = unlist(ret_list), expand.grid(fcst_date = rownames(object$Z_fcst[incl_fcst,,2]),
                                                            time = NA,
                                                            variable = object$names_col,
                                                            quantile = pred_quantiles))
    ret_tidy$time <- rep(nrow(object$Y)+object$n_fcst-max(incl_fcst)+incl_fcst, nrow(ret_tidy)/length(incl_fcst))
    return(ret_tidy)
  }

}
