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
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 100, n_reps = 100)
#' prior_obj <- update_prior(prior_obj, n_fcst = 4)
#'
set_prior <- function(Y, d = NULL, d_fcst = NULL, freq, prior_Pi_AR1 = rep(0, ncol(Y)),
                      lambda1 = 0.2, lambda2 = 1, lambda3 = 10000,
                      prior_psi_mean = NULL, prior_psi_Omega = NULL, n_lags, n_fcst = 0, n_burnin, n_reps,
                      verbose = FALSE, smooth_state = FALSE, check_roots = TRUE) {
  if (!is.matrix(Y)) {
    if (!is.data.frame(Y)) {
      stop(paste0("Y is of class ", class(Y), ", but must be matrix or data frame."))
    } else {
      Y <- as.matrix(Y)
    }
  }

  if (nrow(Y) > max(unlist(apply(Y, 2, function(x) Position(is.na, x, nomatch = nrow(Y)))))) {
    stop("Y: remove final rows containing only NAs.")
  }

  if (min(unlist(apply(Y[, freq == "m", drop = FALSE], 2, function(x) Position(is.na, x, nomatch = 9999999999)))) == 1) {
    stop("Y: monthly variables are NA at the beginning of the sample.")
  }

  intercept_flag <- FALSE
  if (hasArg(d)) {
    if (is.vector(d) && length(d) == 1 && d == "intercept") {
      intercept_flag <- TRUE
      d <- matrix(1, nrow = nrow(Y), 1)
    }
    if (!is.matrix(d)) {
      stop("d is of class", class(d), ", but must be a matrix.")
    }

    if (nrow(Y) != nrow(d)) {
      stop("Y has ", nrow(Y), "rows and d ", nrow(d), "rows, but they must be equal.")
    }

    if (!intercept_flag) {
      if (hasArg(d_fcst)) {
        if (!is.matrix(d_fcst)) {
          stop("d is of class", class(d), ", but must be a matrix.")
        }

        if (ncol(d) != ncol(d_fcst)) {
          stop("d has", ncol(d), " columns and d_fcst", ncol(d_fcst), ".")
        }
      }
    }

  }

  if (hasArg(freq)) {
    if (!is.vector(freq)) {
      stop("freq is of class ", class(freq), ", but it must be a character vector.")
    } else if (!all(freq %in% c("m", "q"))) {
      stop("Elements of freq must be 'm' or 'q'.")
    } else if (length(freq) != ncol(Y)) {
      stop("The length of freq is ", length(freq), ", but Y has ", ncol(Y), " columns.")
    }
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
    } else {
      if (intercept_flag) {
        d_fcst <- matrix(1, nrow = n_fcst, 1)
      }
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
              prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega, n_lags = n_lags, n_fcst = n_fcst, n_burnin = n_burnin,
              n_reps = n_reps, verbose = verbose, smooth_state = smooth_state, check_roots = check_roots, intercept_flag = intercept_flag)
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

#' @rdname set_prior
#'
#' @param prior_obj an object of class \code{mfbvar_prior}
#' @param ... named arguments for prior attributes to update
update_prior <- function(prior_obj, ...)
{
  if(class(prior_obj) != "mfbvar_prior") {
    stop("The object to be updated is not of class mfbvar_prior.")
  }

  tmp <- list(...)
  nam <- names(tmp)
  for(iter in 1:length(tmp)) eval(parse(text=paste0("prior_obj$",nam[iter]," = tmp[[iter]]")))

  if ("d" %in% nam) {
    if (!is.vector(tmp$d) || length(tmp$d) > 1 || tmp$d != "intercept") {
      prior_obj$intercept_flag <- FALSE
    } else {
      prior_obj$intercept_flag <- TRUE
      prior_obj$d <- matrix(1, nrow(prior_obj$Y), 1)
    }
  }

  if ("Y" %in% nam) {
    if (prior_obj$intercept_flag) {
      prior_obj$d <- matrix(1, nrow(prior_obj$Y), 1)
    }
  }

  if (prior_obj$intercept_flag && prior_obj$n_fcst > 0) {
    prior_obj$d_fcst <- matrix(1, prior_obj$n_fcst, 1)
  }

  return(prior_obj)
}

#' Mixed-frequency Bayesian VAR
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @param mfbvar_prior a \code{mfbvar_prior} object
#' @param prior_type either \code{"ss"} (steady-state prior) or \code{"minn"} (Minnesota prior)
#' @param ... additional arguments to \code{update_prior} (if \code{mfbvar_prior} is \code{NULL}, the arguments are passed on to \code{set_prior})
#' @details No (other than very primitive) checks are made on the additional arguments. Thus, for more informative error messages and a safer procedure, create the prior first using \code{set_prior}.
#' @return an object of class \code{mfbvar} and \code{mfbvar_ss}/\code{mfbvar_minn} containing posterior quantities as well as the prior object
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#'
estimate_mfbvar <- function(mfbvar_prior = NULL, prior_type, ...) {
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
    if (!(!is.null(mfbvar_prior$Y) && !is.null(mfbvar_prior$d) && !is.null(mfbvar_prior$prior_psi_mean) && !is.null(mfbvar_prior$prior_psi_Omega) && !is.null(mfbvar_prior$n_lags) && !is.null(mfbvar_prior$n_burnin) && !is.null(mfbvar_prior$n_reps))) {
      test_all <- sapply(mfbvar_prior, is.null)
      test_sub <- test_all[c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags", "n_burnin", "n_reps")]
      stop("Missing elements:", names(test_sub)[which(test_sub)])
    }
    if (mfbvar_prior$n_fcst > 0 && nrow(mfbvar_prior$d_fcst) != mfbvar_prior$n_fcst) {
      stop("d_fcst has ", nrow(mfbvar_prior$d_fcst), " rows, but n_fcst is ", mfbvar_prior$n_fcst, ".")
    }

    priors <- prior_Pi_Sigma(lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1, Y = mfbvar_prior$Y,
                             n_lags = mfbvar_prior$n_lags, prior_nu = n_vars + 2)
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
    if (!(!is.null(mfbvar_prior$Y) && !is.null(mfbvar_prior$n_lags) && !is.null(mfbvar_prior$n_burnin) && !is.null(mfbvar_prior$n_reps))) {
      test_all <- sapply(mfbvar_prior, is.null)
      test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
      stop("Missing elements:", names(test_sub)[which(test_sub)])
    }

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

  if (mfbvar_prior$verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", mfbvar_prior$n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }
  if (prior_type == "ss") {
    burn_in <-  gibbs_sampler(Y = mfbvar_prior$Y, d = mfbvar_prior$d, d_fcst = NULL, freq = mfbvar_prior$freq, prior_Pi_mean = prior_Pi_mean,
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega,
                              n_fcst = 0, n_reps = mfbvar_prior$n_burnin, smooth_state = FALSE, check_roots = mfbvar_prior$check_roots, verbose = mfbvar_prior$verbose)
  } else {
    burn_in <-  gibbs_sampler_minn(Y = mfbvar_prior$Y, freq = mfbvar_prior$freq, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1,
                                   lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2,
                                   lambda3 = mfbvar_prior$lambda3, n_lags = mfbvar_prior$n_lags, n_fcst = 0,
                                   n_reps = mfbvar_prior$n_burnin, smooth_state = FALSE, check_roots = FALSE,
                                   verbose = mfbvar_prior$verbose)
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
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega, n_fcst = mfbvar_prior$n_fcst,
                              n_reps = mfbvar_prior$n_reps+1, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                              init_psi = burn_in$psi[dim(burn_in$psi)[1],], init_Z = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = mfbvar_prior$smooth_state,
                              check_roots = mfbvar_prior$check_roots, verbose = mfbvar_prior$verbose)
  } else {
    main_run <- gibbs_sampler_minn(Y = mfbvar_prior$Y, freq = mfbvar_prior$freq, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1,
                                   lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2, lambda3 = mfbvar_prior$lambda3,
                                   n_lags = mfbvar_prior$n_lags, n_fcst = mfbvar_prior$n_fcst, n_reps = mfbvar_prior$n_reps+1, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                                   init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], init_Z = burn_in$Z[,,dim(burn_in$Z)[3]],
                                   smooth_state = mfbvar_prior$smooth_state, check_roots = mfbvar_prior$check_roots, mfbvar_prior$verbose)
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
                        pred_level = c(0.10, 0.90), ...){
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

  if (!is.null(x$psi)) {

    ss_lower  <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
    ss_median <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
    ss_upper  <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))

    ss <- data.frame(expand.grid(time = plot_range, names_col = names_col), lower = c(ss_lower), median = c(ss_median),
                     upper = c(ss_upper))
    ss$value <- c(as.matrix(x$Y[plot_range,]))
    ss <- na.omit(ss)

    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "grey90"), alpha =1) +
      geom_line(data = ss, aes(y = value))
  }
  if (!is.null(x$n_fcst)) {
    preds <- predict(x, pred_quantiles = c(pred_level[1], 0.5, pred_level[2]))
    fcst <- data.frame(expand.grid(time = (x$n_T+1):(x$n_T+x$n_fcst), names_col = names_col),
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
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "#bdbdbd"), alpha = 1) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = "Note: The forecasts are for the underlying variable.")
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }



  p <- p + facet_wrap(~names_col, scales = "free_y") +
    scale_fill_manual(values = c("#bdbdbd", "grey90"),
                      label = c(paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"),
                                paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
    labs(fill = "Intervals",
         title = "Forecasts and steady state intervals",
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

#' Plotting method for class \code{mfbvar_minn}
#'
#' Method for plotting mfbvar objects.
#' @param x object of class \code{mfbvar_minn}
#' @param plot_start Time period (number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param pred_level A vector with the lower and upper quantiles for the forecast intervals.
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20, n_fcst = 4)
#' mod_minn <- estimate_mfbvar(prior_obj, prior_type = "minn")
#' plot(mod_minn)

plot.mfbvar_minn <- function(x, plot_start = NULL, pred_level = c(0.10, 0.90), ...){
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
    fcst <- data.frame(expand.grid(time = (x$n_T+1):(x$n_T+x$n_fcst), names_col = names_col),
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
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "#bdbdbd"), alpha = 1) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = "Note: The forecasts are for the underlying variable.")
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }



  p <- p + facet_wrap(~names_col, scales = "free_y") +
    scale_fill_manual(values = c("#bdbdbd", "grey90"),
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
  if (is.null(object$n_fcst)) {
    stop("No forecasts exist in the provided object.")
  }

  if (tidy == FALSE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[-(1:object$n_lags),,-1], c(1, 2), quantile, prob = xx))
    names(ret_list) <- paste0("quantile_", pred_quantiles*100)
    return(ret_list)
  } else if (tidy == TRUE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[-(1:object$n_lags),,-1], c(1, 2), quantile, prob = xx))
    ret_tidy <- cbind(value = unlist(ret_list), expand.grid(fcst_date = rownames(object$Z_fcst[-(1:object$n_lags),,2]),
                                                            variable = object$names_col,
                                                            quantile = pred_quantiles))
    return(ret_tidy)
  }

}
