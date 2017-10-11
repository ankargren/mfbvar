#' Set priors for an mfbvar model
#'
#' Create an object storing all information needed for estimation, including data as well as model and prior specifications for both a Minnesota or steady-state prior.
#'
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar d_fcst TRUE
#' @templateVar freq TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar lambda3 TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @templateVar verbose TRUE
#' @templateVar smooth_state TRUE
#' @templateVar check_roots TRUE
#' @template man_template
set_prior <- function(Y, d = NULL, d_fcst = NULL, freq = c("m", "q"), prior_Pi_AR1 = rep(0, ncol(Y)),
                      lambda1 = 0.2, lambda2 = 1, lambda3 = 10000, prior_nu = ncol(Y) + 2,
                      prior_psi_mean = NULL, prior_psi_Omega = NULL, n_lags, n_fcst = 0, n_burnin, n_reps,
                      verbose = FALSE, smooth_state = FALSE, check_roots = TRUE) {
  if (!is.matrix(Y)) {
    if (!is.data.frame(Y)) {
      stop(paste0("Y is of class ", class(Y), ", but must be matrix or data frame."))
    } else {
      Y <- as.matrix(Y)
    }
  }

  if (hasArg(d)) {
    if (!is.matrix(d)) {
      stop("d is of class", class(d), ", but must be a matrix.")
    }

    if (nrow(Y) != nrow(d)) {
      stop("Y has ", nrow(Y), "rows and d ", nrow(d), "rows, but they must be equal.")
    }

    if (hasArg(d_fcst)) {
      if (!is.matrix(d_fcst)) {
        stop("d is of class", class(d), ", but must be a matrix.")
      }

      if (ncol(d) != ncol(d_fcst)) {
        stop("d has", ncol(d), " columns and d_fcst", ncol(d_fcst), ".")
      }
    }
  }

  if (!is.vector(freq)) {
    stop("freq is of class ", class(freq), ", but it must be a character vector.")
  } else if (!all(freq %in% c("m", "q"))) {
    stop("Elements of freq must be 'm' or 'q'.")
  } else if (length(freq) != ncol(Y)) {
    stop("The length of freq is ", length(freq), ", but Y has ", ncol(Y), " columns.")
  }

  if (hasArg(prior_Pi_AR1)) {
    if (is.vector(prior_Pi_AR1)) {
      if (length(prior_Pi_AR1) == 1) {
        warning("prior_Pi_AR1: Recycling ", prior_Pi_AR1, " to use as prior mean for AR(1) coefficients for all variables.\n", call. = FALSE)
        prior_Pi_AR1 <- rep(prior_Pi_AR1, ncol(Y))
      }
    } else {
      stop("prior_Pi_AR1 must be a vector, but is now of class ", class(prior_Pi_AR1), ".")
    }
  } else {
    warning("prior_Pi_AR1: Using 0 as prior mean for AR(1) coefficients for all variables.\n", call. = FALSE)
    prior_Pi_AR1 <- rep(0, ncol(Y))
  }

  if (hasArg(lambda1)) {
    if (!is.vector(lambda1) || length(lambda1) > 1) {
      stop("lambda1 must be a vector with a single element.")
    }
  } else {
    warning("lambda1: Using the default ", lambda1, " as the value for the overall tightness hyperparameter.\n", call. = FALSE)
  }

  if (hasArg(lambda2)) {
    if (!is.vector(lambda2) || length(lambda2) > 1) {
      stop("lambda2 must be a vector with a single element.")
    }
  } else {
    warning("lambda2: Using the default ", lambda2, " as the value for the lag decay hyperparameter.\n", call. = FALSE)
  }

  if (hasArg(lambda3)) {
    if (!is.vector(lambda3) || length(lambda3) > 1) {
      stop("lambda3 must be a vector with a single element.")
    }
  } else {
    warning("lambda3: Using the default ", lambda3, " as the constant's prior variance.\n", call. = FALSE)
  }

  if (hasArg(prior_nu)) {
    if (!is.vector(prior_nu) || length(prior_nu) > 1) {
      stop("prior_nu must be a vector with a single element.")
    }
  } else {
    warning(paste0("prior_nu: Using the default n_vars + 2 = ", prior_nu, " prior for prior_nu.\n"), call. = FALSE)
  }


  if (hasArg(prior_psi_mean)) {
    if (!(is.vector(prior_psi_mean) || is.matrix(prior_psi_mean))) {
      stop("prior_psi_mean must be a vector or matrix with one row or column.")
    }
    if (is.vector(prior_psi_mean)) {
      if (length(prior_psi_mean) != ncol(Y)) {
        stop("prior_psi_mean has ", length(prior_psi_mean), " elements, but there are ", ncol(Y), " variables in Y.")
      }
    }
    if (is.matrix(prior_psi_mean)) {
      if (!any(dim(prior_psi_mean) == 1)) {
        stop("prior_psi_mean must be a vector or matrix with one row or column.")
      } else {
        prior_psi_mean <- c(prior_psi_mean)
      }
    }
  }

  if (hasArg(prior_psi_Omega)) {
    if (!is.matrix(prior_psi_Omega)) {
      stop("prior_psi_Omega must be a matrix.")
    } else {
      if (dim(prior_psi_Omega)[1] != dim(prior_psi_Omega)[2]) {
        stop("prior_psi_Omega must be a positive-definite symmetric matrix.")
      }
    }
  }

  if (hasArg(n_lags)) {
    if (!is.vector(n_lags) || length(n_lags) > 1) {
      stop("n_lags must be a vector with a single element.")
    }
  } else {
    stop("n_lags: No lag length specified.\n")
  }

  if (hasArg(n_fcst)) {
    if (!is.vector(n_fcst) || length(n_fcst) > 1) {
      stop("n_fcst must be a vector with a single element.")
    }
  } else {
    warning("n_fcst: Using the default ", n_fcst, " for the number of forecasts to compute.\n", call. = FALSE)
  }

  if (hasArg(n_burnin)) {
    if (!is.vector(n_burnin) || length(n_burnin) > 1) {
      stop("n_burnin must be a vector with a single element.")
    }
  } else {
    stop("n_burnin: Number of burn-in draws to use not specified.\n")
  }

  if (hasArg(n_reps)) {
    if (!is.vector(n_reps) || length(n_reps) > 1) {
      stop("n_reps must be a vector with a single element.")
    }
  } else {
    stop("n_reps: Number of draws to use in main chain not specified.\n")
  }

  if (!is.logical(smooth_state)) {
    stop("smooth_state: must be logical.\n")
  }

  if (!is.logical(check_roots)) {
    stop("check_roots: must be logical.\n")
  }

  ret <- list(Y = Y, d = d, d_fcst = d_fcst, freq = freq, prior_Pi_AR1 = prior_Pi_AR1, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3,
              prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega, n_lags = n_lags, n_fcst = n_fcst, n_burnin = n_burnin,
              n_reps = n_reps, verbose = verbose, smooth_state = smooth_state, check_roots = check_roots)
  class(ret) <- "mfbvar_prior"

  return(ret)
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
  cat("  d:", ifelse(is.null(object$d), "<missing>", paste0(ncol(object$d), "deterministic variables")),"\n")
  cat("  d_fcst:", ifelse(is.null(object$d_fcst), "<missing>", paste0(nrow(object$d_fcst), "forecasts, ", ncol(object$d), "deterministic variables")),"\n")
  cat("  prior_nu:", ifelse(is.null(object$prior_nu), "<missing> (will rely on default)", object$prior_nu), "\n")
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
#' Update priors for an mfbvar model
#'
#' @param prior_obj an object of class \code{mfbvar\_prior}
#' @param ... named arguments for prior attributes to update
#' @seealso \code{\link{set_prior}}
update_prior <- function(prior_obj, ...)
{
  if(class(prior_obj) != "mfbvar_prior") {
    stop("The object to be updated is not of class mfbvar_prior.")
  }

  tmp <- list(...)
  nam <- names(tmp)
  for(iter in 1:length(tmp)) eval(parse(text=paste0("prior_obj$",nam[iter]," = tmp[[iter]]")))

  return(prior_obj)
}

#' Mixed-frequency Bayesian VAR
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @param mfbvar_prior Test
#' @param prior_type Test 2
#' @param ... test 3

estimate_mfbvar <- function(mfbvar_prior = NULL, prior_type = c("ss", "minn"), ...) {
  if (hasArg(mfbvar_prior)) {
    if (class(mfbvar_prior) != "mfbvar_prior") {
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

  n_vars <- ncol(mfbvar_prior$Y)

  if (prior_type == "ss") {
    if (mfbvar_prior$n_fcst > 0 && nrow(mfbvar_prior$d_fcst) != mfbvar_prior$n_fcst) {
      stop("d_fcst has ", nrow(mfbvar_prior$d_fcst), " rows, but n_fcst is ", mfbvar_prior$n_fcst, ".")
    }

    priors <- prior_Pi_Sigma(lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1, Y = mfbvar_prior$Y,
                             n_lags = mfbvar_prior$n_lags, prior_nu = mfbvar_prior$prior_nu)
    prior_Pi_mean <- priors$prior_Pi_mean
    prior_Pi_Omega <- priors$prior_Pi_Omega
    prior_S <- priors$prior_S

    if (mfbvar_prior$n_fcst > 0) {
      if (is.null(rownames(mfbvar_prior$d_fcst))) {
        names_fcst <- paste0("fcst_", 1:mfbvar_prior$n_fcst)
      } else {
        names_fcst <- rownames(mfbvar_prior$d_fcst)
      }
    } else {
      names_fcst <- NULL
    }

    if (is.null(colnames(mfbvar_prior$d))) {
      names_determ <- paste0("d", 1:ncol(mfbvar_prior$d))
    } else {
      names_determ <- colnames(mfbvar_prior$d)
    }

  } else {
    if (mfbvar_prior$n_fcst > 0) {
      names_fcst <- paste0("fcst_", 1:mfbvar_prior$n_fcst)
    } else {
      names_fcst <- NULL
    }
  }

  if (is.null(rownames(mfbvar_prior$Y))) {
    names_row <- 1:nrow(mfbvar_prior$Y)
  } else {
    names_row <- rownames(mfbvar_prior$Y)
  }

  if (is.null(colnames(mfbvar_prior$Y))) {
    names_col <- 1:col(mfbvar_prior$Y)
  } else {
    names_col <- colnames(mfbvar_prior$Y)
  }


  Lambda <- build_Lambda(mfbvar_prior$freq, mfbvar_prior$n_lags)

  if (mfbvar_prior$verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", mfbvar_prior$n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }
  if (prior_type == "ss") {
    burn_in <-  gibbs_sampler(Y = mfbvar_prior$Y, d = mfbvar_prior$d, d_fcst = NULL, freq = mfbvar_prior$freq, prior_Pi_mean = prior_Pi_mean,
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S, prior_nu = mfbvar_prior$prior_nu,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega,
                              n_fcst = 0, n_reps = mfbvar_prior$n_burnin, smooth_state = FALSE, check_roots = mfbvar_prior$check_roots, verbose = mfbvar_prior$verbose)
  } else {
    burn_in <-  gibbs_sampler_minn(Y = mfbvar_prior$Y, freq = mfbvar_prior$freq, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1,
                                     lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2, lambda3 = mfbvar_prior$lambda3,
                                     n_lags = mfbvar_prior$n_lags, n_fcst = 0, n_reps = mfbvar_prior$n_burnin,
                                     smooth_state = FALSE, check_roots = FALSE, verbose = mfbvar_prior$verbose)
  }

  if (mfbvar_prior$verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", mfbvar_prior$n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\n   Moving on to the main chain with ",
               mfbvar_prior$n_reps, " draws \n\n", ifelse(mfbvar_prior$n_fcst > 0, paste0("   Making forecasts ", mfbvar_prior$n_fcst, " steps ahead"), ""), "\n\n"))
  }

  if (prior_type == "ss") {
    main_run <- gibbs_sampler(Y = mfbvar_prior$Y, d = mfbvar_prior$d, d_fcst = mfbvar_prior$d_fcst, freq = mfbvar_prior$freq, prior_Pi_mean = prior_Pi_mean,
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S, prior_nu = mfbvar_prior$prior_nu,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega, n_fcst = mfbvar_prior$n_fcst,
                              n_reps = mfbvar_prior$n_reps+1, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                              init_psi = burn_in$psi[dim(burn_in$psi)[1],], init_Z = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = mfbvar_prior$smooth_state,
                              check_roots = mfbvar_prior$check_roots, verbose = mfbvar_prior$verbose)
  } else {
    main_run <- gibbs_sampler_minn(Y = mfbvar_prior$Y, freq = mfbvar_prior$freq, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1, lambda1 = mfbvar_prior$lambda1,
                                     lambda2 = mfbvar_prior$lambda2, lambda3 = mfbvar_prior$lambda3, n_lags = mfbvar_prior$n_lags, n_fcst = mfbvar_prior$n_fcst,
                                     n_reps = mfbvar_prior$n_reps+1, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                                     init_Z = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = mfbvar_prior$smooth_state, check_roots = mfbvar_prior$check_roots, mfbvar_prior$verbose)
  }


  if (mfbvar_prior$verbose) {
    time_diff <- Sys.time() - start_burnin
    cat(paste0("\n   Total time elapsed: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
  }

  if (mfbvar_prior$n_fcst > 0) {
    main_run$Z_fcst <- main_run$Z_fcst[,,-1]
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(main_run$n_T-main_run$n_lags+1):main_run$n_T]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col
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

