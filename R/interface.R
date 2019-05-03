#' Set priors for an mfbvar model
#'
#' Create an object storing all information needed for estimation, including data as well as model and prior specifications for both a Minnesota or steady-state prior.
#'
#' @templateVar Y TRUE
#' @templateVar freq TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @param (Only if \code{variance != "iw"}) The cross-variable tightness
#' @templateVar lambda3 TRUE
#' @param lambda4 (Minnesota only) Prior variance of the intercept.
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @param n_thin Store every \code{n_thin}th draw
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @param d (Steady state only) Either a matrix with same number of rows as \code{Y} and \code{n_determ} number of columns containing the deterministic terms or a string \code{"intercept"} for requesting an intercept as the only deterministic
#' term.
#' @param d_fcst (Steady state only) The deterministic terms for the forecasting period (not used if \code{d = "intercept"}).
#' @param prior_psi_mean (Steady state only) Vector of length \code{n_determ*n_vars} with the prior means of the steady-state parameters.
#' @param prior_psi_Omega (Steady state only) Matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)} with the prior covariance of the steady-state parameters.
#' @param prior_phi (Only used with common stochastic volatility) Vector with two elements \code{c(mean, variance)} for the AR(1) parameter in the log-volatility regression
#' @param prior_sigma2 (Only used with common stochastic volatility) Vector with two elements \code{c(mean, df)} for the innovation variance of the log-volatility regression
#' @param n_fac (Only used with factor stochastic volatility) Number of factors to use for the factor stochastic volatility model
#' @param cl (Only used with factor stochastic volatility) Cluster object to use for drawing regression parameters in parallel
#' @param ... (Only used with factor stochastic volatility) Arguments to pass along to \code{\link[factorstochvol]{fsvsample}}. See details.
#' @templateVar verbose TRUE
#' @templateVar check_roots TRUE
#' @template man_template
#' @details The first arguments (\code{Y} through \code{n_reps}) must be set for the model to be estimated irrespective of the choice
#' of prior, but some have default values.
#'
#' For the Minnesota prior, \code{lambda4} must also be set, but it too has a default that it relies on if not specified.
#'
#' For the steady-state prior, the deterministic matrix needs to be supplied, or a string indicating that the intercept should be
#' the only deterministic term. If the latter, also \code{d_fcst} is set to be intercept only. Otherwise, if forecasts are requested
#' (\code{n_fcst > 0}) also \code{d_fcst} needs to be provided. Finally, the prior moments for the steady-state parameters must also be
#' provided.
#'
#' For modeling stochastic volatility by the factor stochastic volatility model, the number of factors to use must be supplied. Further arguments can be passed along to \code{\link[factorstochvol]{fsvsample}}. If arguments are not given, the defaults used are as follows (see \code{\link[factorstochvol]{fsvsample}} for descriptions):
#' \itemize{
#'   \item{\code{priormu}}{\code{ = c(0, 10)}}
#'   \item{\code{priorphiidi}}{\code{ = c(10, 3)}}
#'   \item{\code{priorphifac}}{\code{ = c(10, 3)}}
#'   \item{\code{priorsigmaidi}}{\code{ = 1}}
#'   \item{\code{priorsigmafac}}{\code{ = 1}}
#'   \item{\code{priorfacload}}{\code{ = 1}}
#'   \item{\code{priorng}}{\code{ = c(1, 1)}}
#'   \item{\code{columnwise}}{\code{ = FALSE}}
#'   \item{\code{restrict}}{\code{ = "none"}}
#'   \item{\code{heteroskedastic}}{\code{ = TRUE}}
#'   \item{\code{priorhomoskedastic}}{\code{ = NA}}
#' }
#'
#' The steady-state prior involves inverting the lag polynomial. For this reason, draws in which the largest eigenvalue
#' (in absolute value) of the lag polynomial is greater than 1 are discarded and new draws are made. The maximum number of
#' attempts is 1,000. The components in the output named \code{roots} and \code{num_tries} contain the largest roots and the
#' number of attempts, respectively, if \code{check_roots = TRUE} (the default).
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 100, n_reps = 100)
#' prior_obj <- update_prior(prior_obj, n_fcst = 4)
#' @seealso \code{\link{interval_to_moments}}, \code{\link{print.mfbvar_prior}}, \code{\link{summary.mfbvar_prior}}, \code{\link{estimate_mfbvar}}, \code{\link[factorstochvol]{fsvsample}}
set_prior <- function(Y, freq, prior_Pi_AR1 = rep(0, ncol(Y)), lambda1 = 0.2,
                      lambda2 = 0.5, lambda3 = 1, lambda4 = 10000, n_lags,
                      n_fcst = 0, n_thin = 1, n_burnin, n_reps, d = NULL, d_fcst = NULL,
                      prior_psi_mean = NULL, prior_psi_Omega = NULL, prior_phi = c(0.9, 0.1),
                      prior_sigma2 = c(0.01, 4), n_fac = NULL,
                      cl = NULL, verbose = FALSE, check_roots = FALSE, ...) {
  prior_call <- mget(names(formals())[names(formals()) != "..."], sys.frame(sys.nframe()))
  prior_call$supplied_args <- names(as.list(match.call()))[-1]
  ellipsis <- list(...)
  fsv_names <- names(ellipsis)
  fsv_arguments <- c("priormu", "priorphiidi", "priorphifac", "priorsigmaidi", "priorsigmafac", "priorfacload", "priorng", "columnwise", "restrict", "heteroskedastic", "priorhomoskedastic")
  if (!all(fsv_names %in% fsv_arguments)) {
    unused_names <- setdiff(fsv_names, fsv_arguments)
    warning(sprintf("The following arguments passed along to fsvsample are unused: %s", ifelse(unused_names == "", "[unnamed component]", unused_names)))
  }
  prior_call <- append(prior_call, ellipsis[fsv_names %in% fsv_arguments])
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
      col_class <- sapply(prior_obj$Y, class)
      if (all(col_class == "numeric")) {
        prior_obj$Y <- as.matrix(prior_obj$Y)
        if (inherits(prior_obj$Y, "tbl")) {
          rownames(prior_obj$Y) <- 1:nrow(prior_obj$Y)
        } else {
          if (is.null(rownames(prior_obj$Y))) {
            rownames(prior_obj$Y) <- 1:nrow(prior_obj$Y)
          }
        }
      } else if (sum(!(col_class == "numeric")) == 1) {
          date_tmp <- prior_obj$Y[[which(col_class != "numeric")]]
          date_tmp <- as.character(as.Date(date_tmp))
          prior_obj$Y <- as.matrix(prior_obj$Y[, which(col_class == "numeric")])
          if (!any(is.na(date_tmp))) {
            rownames(prior_obj$Y) <- date_tmp
          }
      }
      else {
        stop(sprintf("The data frame contains %d non-numeric columns. Please include at most one non-numeric column that can be coerced to dates.", sum(!(col_class == "numeric"))))
      }
    }
  } else {
    if (is.null(rownames(prior_obj$Y))) {
      rownames(prior_obj$Y) <- 1:nrow(prior_obj$Y)
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
        prior_obj$prior_Pi_AR1 <- rep(prior_obj$prior_Pi_AR1, ncol(prior_obj$Y))
      }
    } else {
      stop("prior_Pi_AR1 must be a vector, but is now of class ", class(prior_obj$prior_Pi_AR1), ".")
    }
  } else {
    prior_obj$prior_Pi_AR1 <- rep(0, ncol(prior_obj$Y))
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "prior_Pi_AR1")
  }


  if ("lambda1" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda1) || length(prior_obj$lambda1) > 1) {
      stop("lambda1 must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda1")
  }

  if ("lambda2" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda2) || length(prior_obj$lambda2) > 1) {
      stop("lambda2 must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda2")
  }

  if ("lambda3" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda3) || length(prior_obj$lambda3) > 1) {
      stop("lambda3 must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda3")
  }

  if ("lambda4" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$lambda4) || length(prior_obj$lambda4) > 1) {
      stop("lambda4 must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "lambda4")
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
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "n_fcst")
  }

  if ("n_thin" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_thin) || length(prior_obj$n_thin) > 1) {
      stop("n_thin must be a vector with a single element.")
    }
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

  if (!is.logical(prior_obj$check_roots)) {
    stop("check_roots: must be logical.\n")
  }

  if ("prior_phi" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$prior_phi) || length(prior_obj$prior_phi) != 2) {
      stop("prior_phi must be a vector with two numeric elements.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "prior_phi")
  }

  if ("prior_sigma2" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$prior_sigma2) || length(prior_obj$prior_sigma2) != 2) {
      stop("prior_sigma2 must be a vector with two numeric elements.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "prior_sigma2")
  }


  if (!is.null(prior_obj$n_fac)) {
    if (!is.numeric(prior_obj$n_fac) | !is.atomic(prior_obj$n_fac)) {
      stop("The number of factors is not a numeric scalar value.")
    }

    if (!inherits(prior_obj$cl, "cluster") && !is.null(prior_obj$cl)) {
      stop(sprintf("cl should be a cluster object, but is %s", class(prior_obj$cl)))
    }

    if ("priormu" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priormu) && length(prior_obj$priormu) == 2)) {
        stop(sprintf("priormu should be a numeric vector of length 2, but is %s of length %d", class(prior_obj$priormu), length(prior_obj$priormu)))
      }
    } else {
      prior_obj$priormu <- c(0, 10)
    }

    if ("priorphiidi" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorphiidi) && length(prior_obj$priorphiidi) == 2)) {
        stop(sprintf("priorphiidi should be a numeric vector of length 2, but is %s of length %d", class(prior_obj$priorphiidi), length(prior_obj$priorphiidi)))
      }
    } else {
      prior_obj$priorphiidi <- c(10, 3)
    }

    if ("priorphifac" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorphifac) && length(prior_obj$priorphifac) == 2)) {
        stop(sprintf("priorphifac should be a numeric vector of length 2, but is %s of length %d", class(prior_obj$priorphifac), length(prior_obj$priorphifac)))
      }
    } else {
      prior_obj$priorphifac <- c(10, 3)
    }

    if ("priorsigmaidi" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorsigmaidi) && is.atomic(prior_obj$priorsigmaidi) && length(prior_obj$priorsigmaidi) %in% c(1, ncol(prior_obj$Y)))) {
        stop(sprintf("priorsigmaidi should be a numeric vector with 1 or n_vars elements, but is %s with %d elements", class(prior_obj$priorsigmaidi), length(prior_obj$priorsigmaidi)))
      }
    } else {
      prior_obj$priorsigmaidi <- 1
    }

    if ("priorsigmafac" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorsigmafac) && is.atomic(prior_obj$priorsigmafac) && length(prior_obj$priorsigmafac) %in% c(1, ncol(prior_obj$n_fac)))) {
        stop(sprintf("priorsigmafac should be a numeric vector with 1 or n_vars elements, but is %s with %d elements", class(prior_obj$priorsigmafac), length(prior_obj$priorsigmaidi)))
      }
    } else {
      prior_obj$priorsigmafac <- 1
    }

    if ("priorfacload" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$supplied_args$priorfacload) && (length(prior_obj$supplied_args$priorfacload) == 1 || dim(prior_obj$supplied_args$priorfacload) == c(ncol(prior_obj$Y), prior_obj$factors)))) {
        stop(sprintf("priorfacload should be a scalar value or an n_vars x n_fac matrix, but is %s with %d elements", class(prior_obj$priorfacload), length(prior_obj$priorfacload)))
      }
    } else {
      prior_obj$priorfacload <- 1
    }

    if ("priorng" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorng) && length(prior_obj$priorng) == 2)) {
        stop(sprintf("priorng should be a numeric vector of length 2, but is %s of length %d", class(prior_obj$priorng), length(prior_obj$priorng)))
      }
    } else {
      prior_obj$priorng <- c(1, 1)
    }

    if ("columnwise" %in% prior_obj$supplied_args) {
      if (!(is.logical(prior_obj$columnwise) && length(prior_obj$priorng) == 1)) {
        stop(sprintf("columnwise should be a single logical value, but is %s of length %d", class(prior_obj$columnwise), length(prior_obj$columnwise)))
      }
    } else {
      prior_obj$columnwise <- FALSE
    }

    if ("restrict" %in% prior_obj$supplied_args) {
      if (!(is.character(prior_obj$restrict) && length(prior_obj$priorng) == 1)) {
        stop(sprintf("restrict should be a single string, but is %s of length %d", class(prior_obj$restrict), length(prior_obj$restrict)))
      }
    } else {
      prior_obj$restrict <- "none"
    }

    if ("heteroskedastic" %in% prior_obj$supplied_args) {
      if (!(is.logical(prior_obj$heteroskedastic) && length(prior_obj$priorng) %in% c(1, 2, ncol(prior_obj$Y)+prior_obj$n_fac))) {
        stop(sprintf("heteroskedastic should be a vector of 1, 2, or n_vars + n_fac logical values, but is %s of length %d", class(prior_obj$heteroskedastic), length(prior_obj$heteroskedastic)))
      }
    } else {
      prior_obj$heteroskedastic <- TRUE
    }

    if (any(!prior_obj$heteroskedastic)) {
      if ("priorhomoskedastic" %in% prior_obj$supplied_args) {
        if (!(is.numeric(prior_obj$priorhomoskedastic) && is.matrix(prior_obj$priorhomoskedastic) && dim(prior_obj$priorhomoskedastic) == c(ncol(prior_obj$Y)+prior_obj$n_fac, 2))) {
          stop(sprintf("priorhomoskedastic should be a matrix of dimensions (n_vars + n_fac) x 2, but is %s of length %d", class(prior_obj$priorhomoskedastic), length(prior_obj$priorhomoskedastic)))
        }
      } else {
        prior_obj$priorhomoskedastic <- c(1.1, 1.1)
      }
    } else {
      prior_obj$priorhomoskedastic <- c(1.1, 1.1)
    }
  } else if (is.null(prior_obj$n_fac) && any(prior_obj$supplied_args %in% c("priormu", "priorphiidi", "priorphifac", "priorsigmaidi", "priorsigmafac",
                   "priorfacload", "priorng", "columnwise", "restrict", "heteroskedastic", "priorhomoskedastic"))) {
    stop("Please set the number of factors before attempting to pass additional arguments along to fsvsim.")
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
#'   priors can be used with the current specification. This check is minimal in
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
  cat("Checking if steady-state prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$prior_psi_Omega) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n")
  }

  cat("Checking if Minnesota prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
  }

  cat("Checking if common stochastic volatility can be used... ")
  if (length(x$prior_phi) == 2 && length(x$prior_sigma2) == 2) {
    cat("TRUE\n\n")
  } else {
    switch(paste0(as.numeric(is.null(x$prior_phi)), as.numeric(is.null(x$prior_sigma2))),
           "01" = cat("FALSE\n Missing element: prior_sigma2 \n\n"),
           "10" = cat("FALSE\n Missing element: prior_phi \n\n"),
           "00" = cat("FALSE\n Missing elements: prior_phi, prior_sigma2 \n\n"))

  }

  cat("Checking if factor stochastic volatility can be used... ")
  if (!is.null(x$n_fac)) {
    cat("TRUE\n\n")
  } else {
    cat("FALSE\n Missing element: n_fac \n\n")
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
  cat("Main specification:\n")
  cat("  Y:", ncol(object$Y), "variables,", nrow(object$Y), "time points\n")
  cat(sprintf("  freq: %d monthly and %d quarterly %s\n", sum(object$freq == "m"), sum(object$freq == "q"), ifelse(sum(object$freq == "q") == 1, "variable", "variables")))
  if (length(object$prior_Pi_AR1)<=5) {
    disp_prior_Pi_AR1 <- object$prior_Pi_AR1
  } else {
    if (length(unique(object$prior_Pi_AR1)) == 1) {
      disp_prior_Pi_AR1 <- object$prior_Pi_AR1[1]
    } else {
      disp_prior_Pi_AR1 <- sprintf("vector with %d elements", length(object$prior_Pi_AR1))
    }
  }
  cat("  prior_Pi_AR1:", disp_prior_Pi_AR1, "\n")
  cat("  lambda1:", object$lambda1, "\n")
  cat("  lambda2:", object$lambda2, "\n")
  cat("  lambda3:", object$lambda3, "\n")
  cat("  lambda4:", ifelse(is.null(object$lambda4), "<missing> (will rely on default)", object$lambda4), "\n")
  cat("  n_lags:", object$n_lags, "\n")
  cat("  n_fcst:", object$n_fcst, "\n")
  cat("  n_burnin:", object$n_burnin, "\n")
  cat("  n_reps:", object$n_reps, "\n")
  cat("----------------------------\n")
  cat("Steady-state-specific elements:\n")
  cat("  d:", ifelse(is.null(object$d), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(ncol(object$d), "deterministic variables"))),"\n")
  cat("  d_fcst:", ifelse(is.null(object$d_fcst), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(nrow(object$d_fcst), "forecasts, ", ncol(object$d), "deterministic variables"))),"\n")
  cat("  prior_psi_mean:", ifelse(is.null(object$prior_psi_mean), "<missing>", "prior mean vector of steady states"), "\n")
  cat("  prior_psi_Omega:", ifelse(is.null(object$prior_psi_Omega), "<missing>", "prior covariance matrix of steady states"), "\n")
  cat("----------------------------\n")
  cat("Factor stochastic volatility-specific elements:\n")
  cat("  n_fac:", ifelse(is.null(object$n_fac), "<missing>", object$n_fac), "\n")
  cat("  cl:", ifelse(is.null(object$cl), "<missing>", sprintf("%s with %d workers", class(object$cl)[1], length(object$cl))), "\n")
  if ("priormu" %in% object$supplied_args) {
    cat("  priormu:", object$priormu, "\n")
  }
  if ("priorphiidi" %in% object$supplied_args) {
    cat("  priorphiidi:", object$priorphiidi, "\n")
  }
  if ("priorphifac" %in% object$supplied_args) {
    cat("  priorphifac:", object$priorphifac, "\n")
  }
  if ("priorsigmaidi" %in% object$supplied_args) {
    if (length(object$priorsigmaidi) == 1) {
      cat("  priorsigmaidi:", object$priorsigmaidi, "\n")
    } else {
      cat("  priorsigmaidi: vector with", length(object$priorsigmaidi), "elements \n")
    }
  }
  if ("priorsigmafac" %in% object$supplied_args) {
    if (length(object$priorsigmafac) == 1) {
      cat("  priorsigmafac:", object$priorsigmafac, "\n")
    } else {
      cat("  priorsigmafac: vector with", length(object$priorsigmafac), "elements \n")
    }
  }
  if ("priorfacload" %in% object$supplied_args) {
    if (length(object$priorfacload) == 1) {
      cat("  priorfacload:", object$priorfacload, "\n")
    } else {
      cat("  priorfacload:", paste(dim(object$priorfacload), collapse = " x "), "matrix\n")
    }
  }
  if ("priorng" %in% object$supplied_args) {
    cat("  priorng:", object$priorng, "\n")
  }
  if ("columnwise" %in% object$supplied_args) {
    cat("  columnwise:", object$columnwise, "\n")
  }
  if ("restrict" %in% object$supplied_args) {
    cat("  restrict:", object$restrict, "\n")
  }
  if ("heteroskedastic" %in% object$supplied_args) {
    if (length(object$heteroskedastic) <= 2) {
      cat("  heteroskedastic:", object$heteroskedastic, "\n")
    } else {
      cat("  heteroskedastic: vector with", ncol(object$Y)+object$n_fac, "logical values\n")
    }
  }
  if ("priorhomoskedastic" %in% object$supplied_args) {
    cat("  priorhomoskedastic:", paste(dim(object$priorhomoskedastic), collapse = " x "), "matrix\n")
  }
  cat("----------------------------\n")
  cat("Other:\n")
  cat("  verbose:", object$verbose, "\n")
  cat("  check_roots:", object$check_roots, "\n")

}


#' Mixed-frequency Bayesian VAR
#'
#' The main function for estimating a mixed-frequency BVAR.
#'
#' @param mfbvar_prior a \code{mfbvar_prior} object
#' @param prior either \code{"ss"} (steady-state prior) or \code{"minn"} (Minnesota prior)
#' @param variance form of the error variance-covariance matrix: \code{"iw"} for the inverse Wishart prior, \code{"csv"} for common stochastic volatility or \code{"fsv"} for factor stochastic volatility
#' @param ... additional arguments to \code{update_prior} (if \code{mfbvar_prior} is \code{NULL}, the arguments are passed on to \code{set_prior})
#' @return An object of class \code{mfbvar}, \code{mfbvar_<prior>} and \code{mfbvar_<prior>_<variance>} containing posterior quantities as well as the prior object
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{predict.mfbvar}}, \code{\link{plot.mfbvar_minn}},
#' \code{\link{plot.mfbvar_ss}}, \code{\link{varplot}}, \code{\link{summary.mfbvar}}
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden, freq = c(rep("m", 4), "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' @return For all choices of \code{prior} and \code{variance}, the returned object contains:
#' \item{Pi}{Array of dynamic coefficient matrices; \code{Pi[,, r]} is the \code{r}th draw}
#' \item{Z}{Array of monthly processes; \code{Z[,, r]} is the \code{r}th draw}
#' \item{Z_fcst}{Array of monthly forecasts; \code{Z_fcst[,, r]} is the \code{r}th forecast. The first \code{n_lags}
#' rows are taken from the data to offer a bridge between observations and forecasts and for computing nowcasts (i.e. with ragged edges).}
#'
#' If \code{prior = "ss"}, it also includes:
#' \item{psi}{Matrix of steady-state parameter vectors; \code{psi[r,]} is the \code{r}th draw}
#' \item{roots}{The maximum eigenvalue of the lag polynomial (if \code{check_roots = TRUE})}
#' \item{num_tries}{The number of attempts for drawing a stationary \eqn{\Pi} (if \code{check_roots = TRUE})}
#'
#' If \code{variance = "iw"}, it also includes:
#' \item{Sigma}{Array of error covariance matrices; \code{Sigma[,, r]} is the \code{r}th draw}
#'
#' #' If \code{variance = "csv"}, it also includes:
#' \item{Sigma}{Array of error covariance matrices; \code{Sigma[,, r]} is the \code{r}th draw}
#' \item{phi}{Vector of AR(1) parameters for the log-volatility regression; \code{phi[r]} is the \code{r}th draw}
#' \item{sigma}{Vector of error standard deviations for the log-volatility regression; \code{sigma[r]} is the \code{r}th draw}#'
#' \item{f}{Matrix of log-volatilities; \code{f[r, ]} is the \code{r}th draw}
#'
#' If \code{variance = "fsv"}, it also includes:
#' \item{facload}{Array of factor loadings; \code{facload[,, r]} is the \code{r}th draw}
#' \item{latent}{Array of latent log-volatilities; \code{latent[,, r]} is the \code{r}th draw}
#' \item{mu}{Matrix of means of the log-volatilities; \code{mu[, r]} is the \code{r}th draw}
#' \item{phi}{Matrix of AR(1) parameters for the log-volatilities; \code{phi[, r]} is the \code{r}th draw}
#' \item{sigma}{Matrix of innovation variances for the log-volatilities; \code{sigma[, r]} is the \code{r}th draw}
#' @references
#' Schorfheide, F., & Song, D. (2015) Real-Time Forecasting With a Mixed-Frequency VAR. \emph{Journal of Business & Economic Statistics}, 33(3), 366--380. \url{http://dx.doi.org/10.1080/07350015.2014.954707}\cr
#' Ankargren, S., Unosson, M., & Yang, Y. (2018) A Mixed-Frequency Bayesian Vector Autoregression with a Steady-State Prior. Working Paper, Department of Statistics, Uppsala University No. 2018:3.
estimate_mfbvar <- function(mfbvar_prior = NULL, prior, variance = "iw", ...) {
  time_out <- Sys.time()
  args <- list(...)
  if (hasArg(mfbvar_prior)) {
    if (!inherits(mfbvar_prior, "mfbvar_prior")) {
      stop("mfbvar_prior must be of class mfbvar_prior.")
    } else {
      if (length(args) > 0) {
        mfbvar_prior <- update_prior(mfbvar_prior, ...)
      }
    }
  } else {
    mfbvar_prior <- set_prior(...)
  }

  if (hasArg(prior_type)) {
    warning("The argument 'prior_type' is deprecated (starting in 0.5.0). Use 'prior' instead.", call. = FALSE)
    prior <- args$prior_type
  }

  if (!(prior %in% c("ss", "minn"))) {
    stop("prior must be 'ss' or 'minn'.")
  }
  if (!(variance %in% c("iw", "fsv", "csv"))) {
    stop("volatility must be 'iw', 'csv' or 'fsv'.")
  }

  class(mfbvar_prior) <- c(class(mfbvar_prior), sprintf("mfbvar_%s_%s", prior, variance), sprintf("mfbvar_%s", prior), sprintf("mfbvar_%s", variance))

  if (mfbvar_prior$verbose) {
    cat(paste0("##############################################\nRunning the burn-in sampler with ", mfbvar_prior$n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }

  time_out <- c(time_out, Sys.time())
  burn_in <- mcmc_sampler(update_prior(mfbvar_prior, n_fcst = 0), n_reps = mfbvar_prior$n_burnin, n_thin = mfbvar_prior$n_burnin)

  if (mfbvar_prior$verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", mfbvar_prior$n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\nMoving on to the main chain with ",
               mfbvar_prior$n_reps, " draws \n\n", ifelse(mfbvar_prior$n_fcst > 0, paste0("   Making forecasts ", mfbvar_prior$n_fcst, " steps ahead"), ""), "\n\n"))
  }

  time_out <- c(time_out, Sys.time())
  main_run <-  mcmc_sampler(mfbvar_prior, n_reps = mfbvar_prior$n_reps+1, init = burn_in$init)
  time_out <- c(time_out, Sys.time())
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

  dimnames(main_run$Z) <- list(time = names_row[(nrow(mfbvar_prior$Y)-nrow(main_run$Z)+1):nrow(mfbvar_prior$Y)],
                               variable = names_col,
                               iteration = 1:mfbvar_prior$n_reps)

  if (variance == "iw") {
    dimnames(main_run$Sigma) <- list(names_col,
                                     names_col,
                                     iteration = 1:mfbvar_prior$n_reps)
  }


  if (prior == "ss") {
    if (is.null(colnames(mfbvar_prior$d))) {
      names_determ <- paste0("d", 1:ncol(mfbvar_prior$d))
    } else {
      names_determ <- colnames(mfbvar_prior$d)
    }
    rownames(mfbvar_prior$d) <- rownames(mfbvar_prior$Y)
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

  class(main_run) <- c("mfbvar", sprintf("mfbvar_%s_%s", prior, variance), sprintf("mfbvar_%s", prior))
  time_out <- c(time_out, Sys.time())
  main_run$time_out <- time_out
  main_run$variance <- variance
  main_run$prior <- prior
  return(main_run)
}


#' Printing method for class mfbvar
#'
#' Method for printing \code{mfbvar} objects.
#'
#' @param x object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' mod_minn

print.mfbvar <- function(x, ...){
  ss <- ifelse(x$prior == "ss", "steady-state ", "")
  var_type <- switch(x$variance,
                iw = "Inverse Wishart",
                fsv = sprintf("Factor stochastic volatility (%d factors)", x$mfbvar_prior$n_fac),
                csv = "Common stochastic volatility")
  cat(paste0(sprintf("Mixed-frequency %sBVAR with:\n", ss), ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "), "\nError covariance matrix: ", var_type, "\n",
             x$n_lags, " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecasted\n",
             x$n_reps, " draws used in main chain"))
}

#' Summary method for class mfbvar
#'
#' Method for summarizing \code{mfbvar} objects.
#'
#' @param x object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' summary(mod_minn)

summary.mfbvar <- function(x, ...){
  ss <- ifelse(x$prior == "ss", "steady-state ", "")
  var_type <- switch(x$variance,
                     iw = "Inverse Wishart",
                     fsv = sprintf("Factor stochastic volatility (%d factors)", x$mfbvar_prior$n_fac),
                     csv = "Common stochastic volatility")
  cat(paste0(sprintf("Mixed-frequency %sBVAR with:\n", ss), ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "), "\nError covariance matrix: ", var_type, "\n",
             x$n_lags, " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecasted\n",
             x$n_reps, " draws used in main chain"))
}

#' Plotting methods for posterior mfbvar objects
#'
#' Methods for plotting posterior mfbvar objects (\code{mfbvar_minn} and \code{mfbvar_ss}).
#' @param x object of class \code{mfbvar_minn} or \code{mfbvar_ss}
#' @param fcst_start Date of the first forecast; if dates are available for the data used for obtaining \code{x}, these will be used.
#' @param plot_start Time period (date or number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param pred_bands Single number (between \code{0.0} and \code{1.0}) giving the coverage level of forecast intervals.
#' @param ss_bands (Steady-state prior only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior steady-state intervals.
#' @param var_bands (\code{varplot} only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior intervals for the error standard deviations.
#' @param nrow_facet an integer giving the number of rows to use in the facet
#' @param ... Currently not in use.
#' @name plot-mfbvar
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], d = "intercept",
#'                        freq = c("m", "q"), n_lags = 4, n_burnin = 20, n_reps = 20,
#'                        n_fcst = 4, n_fac = 1)
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
#' mod_ss <- estimate_mfbvar(prior_obj, prior = "ss", variance = "fsv")
#' plot(mod_ss)
#' varplot(mod_ss)

#' @rdname plot-mfbvar
plot.mfbvar_ss <- function(x, fcst_start = NULL, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ss_bands = 0.95, ...){


  if (is.null(fcst_start)) {
    row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
    if (inherits(row_names, "error")) {
      stop("To plot the forecasts, either fcst_start must be supplied or the rownames of Y be dates (YYYY-MM-DD).")
    }
    fcst_start <- as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)
  } else {
    fcst_start <- tryCatch(as.Date(fcst_start), error = function(cond) cond)
    if (inherits(fcst_start, "error")) {
      stop("Unable to convert fcst_start to a date.")
    }
  }

  plot_range_names <- fcst_start %m+% months(-x$n_T:(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(x$n_T-x$n_fcst*5, 0):x$n_T
    }  else {
      plot_range <- 1:x$n_T
    }
  } else {
    plot_start <- tryCatch(as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):x$n_T
    } else {
      stop("Unable to convert plot_start to a date.")
    }
  }


  if (is.null(ss_bands)) {
    ss_level <- c(0.025, 0.975)
  } else {
    ss_level <- c(0.5-ss_bands/2, 0.5+ss_bands/2)
  }
  if (is.null(pred_bands)) {
    pred_level <- c(0.10, 0.90)
  } else {
    pred_level <- c(0.5-pred_bands/2, 0.5+pred_bands/2)
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
    ss <- data.frame(expand.grid(time = plot_range[1]:(plot_range[length(plot_range)]+x$n_fcst), variable = names_col), lower = c(ss_lower), median = c(ss_median),
                     upper = c(ss_upper))
    ss$value <- c(rbind(as.matrix(x$Y[plot_range,]), matrix(NA, nrow = x$n_fcst, ncol = x$n_vars)))
    ss_excl <- c(rbind(!is.na(as.matrix(x$Y[plot_range,])), matrix(TRUE, nrow = x$n_fcst, ncol = x$n_vars)))
    #ss <- na.omit(ss)

    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "#bdbdbd"),  alpha =1) +
      geom_line(data = na.omit(ss), aes(y = value))
  }
  if (x$n_fcst > 0) {
    preds <- predict(x, aggregate_fcst = aggregate_fcst, fcst_start = fcst_start, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(variable = names_col[i],
                                     time = (1:nrow(x$Y))[last_pos[i]],
                                     fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]))
    }
    fcst <- mutate(fcst, variable = factor(variable, levels = names_col, labels = names_col))
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "grey90"), linetype = "dotted", color = "black",
                         alpha = 0.75) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = ifelse(aggregate_fcst,
                            "Note: The forecasts for the quarterly variables have been aggregated to the quarterly frequency.",
                            "Note: The forecasts are for the underlying variable."))
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }

  if (x$n_fcst > 0) {
    p <- p + scale_fill_manual(values = c("grey90" = "grey90", "#bdbdbd" = "#bdbdbd"),
                               label = c("grey90" = paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"),
                               "#bdbdbd" = paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Forecasts and posterior steady-state intervals",
           y = "Value",
           x = "Time") +
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd", "grey90"),
                                                     linetype = c("blank", "dotted"))))
  } else {
    p <- p + scale_fill_manual(values = c("#bdbdbd"),
                               label = c(paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Steady-state intervals",
           y = "Value",
           x = "Time")+
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd"),
                                                     linetype = c("blank"))))
  }


  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position="bottom")
  breaks <- ggplot_build(p)$layout$coord$labels(ggplot_build(p)$layout$panel_params)[[1]]$x.labels
  if (any(as.numeric(breaks)>plot_range[length(plot_range)])) {
    break_labels <- c(as.character(plot_range_names[as.numeric(breaks)[as.numeric(breaks)<=plot_range[length(plot_range)]]]),
                      as.character(preds$fcst_date[min(which(preds$time == breaks[as.numeric(breaks)>plot_range[length(plot_range)]]))]))
  } else {
    break_labels <- plot_range_names[as.numeric(breaks)]
  }
  p + scale_x_continuous(breaks = as.numeric(breaks),
                         labels = break_labels) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

#' @rdname plot-mfbvar
plot.mfbvar_minn <- function(x, fcst_start = NULL, aggregate_fcst = TRUE, plot_start = NULL,
                             pred_bands = 0.8, nrow_facet = NULL, ...){
  if (is.null(fcst_start)) {
    row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
    if (inherits(row_names, "error")) {
      stop("To plot the forecasts, either fcst_start must be supplied or the rownames of Y be dates (YYYY-MM-DD).")
    }
    fcst_start <- as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)
  } else {
    fcst_start <- tryCatch(as.Date(fcst_start), error = function(cond) cond)
    if (inherits(fcst_start, "error")) {
      stop("Unable to convert fcst_start to a date.")
    }
  }

  plot_range_names <- fcst_start %m+% months(-x$n_T:(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(x$n_T-x$n_fcst*5, 0):x$n_T
    }  else {
      plot_range <- 1:x$n_T
    }
  } else {
    plot_start <- tryCatch(as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):x$n_T
    } else {
      stop("Unable to convert plot_start to a date.")
    }
  }




  if (is.null(pred_bands)) {
    pred_level <- c(0.10, 0.90)
  } else {
    pred_level <- c(0.5-pred_bands/2, 0.5+pred_bands/2)
  }

  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  p <- ggplot(mapping = aes(x = time))


  ss <- data.frame(expand.grid(time = plot_range[1]:(plot_range[length(plot_range)]+x$n_fcst), variable = names_col))
  ss$value <- c(rbind(as.matrix(x$Y[plot_range,]), matrix(NA, nrow = x$n_fcst, ncol = x$n_vars)))
  ss_excl <- c(rbind(!is.na(as.matrix(x$Y[plot_range,])), matrix(TRUE, nrow = x$n_fcst, ncol = x$n_vars)))
  #ss <- na.omit(ss)

  p <- p +
    geom_line(data = na.omit(ss), aes(y = value))

  if (x$n_fcst > 0) {
    preds <- predict(x, aggregate_fcst = aggregate_fcst, fcst_start = fcst_start, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(variable = names_col[i],
                                     time = (1:nrow(x$Y))[last_pos[i]],
                                     fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]))
    }
    fcst <- mutate(fcst, variable = factor(variable, levels = names_col, labels = names_col))
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "grey90"), linetype = "dotted", color = "black",
                         alpha = 0.75) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = ifelse(aggregate_fcst,
                            "Note: The forecasts for the quarterly variables have been aggregated to the quarterly frequency.",
                            "Note: The forecasts are for the underlying variable."))
  }

  if (x$n_fcst > 0) {
    p <- p + scale_fill_manual(values = c("grey90" = "grey90"),
                               label = c("grey90" = paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Forecast intervals",
           y = "Value",
           x = "Time") +
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd"),
                                                     linetype = c("blank"))))
  } else {
    p <- p +
      labs(y = "Value",
           x = "Time")
  }


  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position="bottom")
  breaks <- ggplot_build(p)$layout$coord$labels(ggplot_build(p)$layout$panel_params)[[1]]$x.labels
  if (any(as.numeric(breaks)>plot_range[length(plot_range)])) {
    break_labels <- c(as.character(plot_range_names[as.numeric(breaks)[as.numeric(breaks)<=plot_range[length(plot_range)]]]),
                      as.character(preds$fcst_date[min(which(preds$time == breaks[as.numeric(breaks)>plot_range[length(plot_range)]]))]))
  } else {
    break_labels <- plot_range_names[as.numeric(breaks)]
  }
  p + scale_x_continuous(breaks = as.numeric(breaks),
                         labels = break_labels) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

#' @rdname plot-mfbvar
varplot <- function(x, variables = colnames(x$Y), var_bands = 0.95, nrow_facet = NULL, ...) {
  if (!inherits(x, c("mfbvar_csv", "mfbvar_fsv"))) {
    stop("The fitted model does not have a time-varying error covariance matrix.")
  }
  if (inherits(x, "mfbvar_csv")) {
    sv_type <- "csv"
  }
  if (inherits(x, "mfbvar_fsv")) {
    sv_type <- "fsv"
  }
  n_reps <- dim(x$latent)[3]
  n_T <- nrow(x$latent)
  n_vars <- ncol(x$Y)
  n_plotvars <- length(variables)
  n_lags <- x$n_lags
  variances <- array(0, dim = c(n_T, n_plotvars, n_reps))
  if (is.character(variables)) {
    variables_num <- which(variables == colnames(x$Y))
  } else {
    variables_num <- variables
  }
  if (sv_type == "fsv") {
    n_fac <- x$mfbvar_prior$n_fac
    for (i in 1:n_reps) {
      for (tt in 1:n_T) {
        variances[tt,,i] <- sqrt(diag(matrix(x$facload[variables_num,,i], n_plotvars, n_fac) %*% diag(exp(x$latent[tt, (n_vars+1):(n_vars+n_fac), i]), n_fac) %*% t(matrix(x$facload[variables_num,,i], n_plotvars, n_fac)))+exp(x$latent[tt, variables_num, i]))
      }
    }
  }
  if (sv_type == "csv") {
    for (i in 1:n_T) {
      for (j in 1:n_reps) {
        variances[i,,j] = exp(0.5*f[j,i])*sqrt(diag(Sigma[,,j])[variables_num])
      }
    }
  }

  date <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(date, "error")) {
    if (is.null(rownames(x$Y))) {
      date <- 1:nrow(x$Y)
    } else {
      date <- as.numeric(rownames(x$Y))
    }
  }
  p <- tibble(date = rep(date[(n_lags+1):(n_lags+n_T)], n_plotvars),
         lower = c(apply(variances, 1:2, quantile, prob = (1-var_bands)/2)),
         mean = c(apply(variances, 1:2, mean)),
         upper = c(apply(variances, 1:2, quantile, prob = 1-(1-var_bands)/2)),
         variable = rep(variables, each = n_T)) %>%
    ggplot(aes(x = date, y = mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90", linetype = "dotted",
                color = "black", alpha = 0.75) +
    geom_line() +
    theme_minimal() +
    labs(y = "Error standard deviation",
         x = "Time")
  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }
  p
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
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' predict(mod_minn)
#' predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE)
predict.mfbvar <- function(object, fcst_start = NULL, aggregate_fcst = TRUE, pred_bands = 0.8, ...) {

  if (object$n_fcst==0) {
    stop("No forecasts exist in the provided object.")
  }
  if (!is.null(fcst_start)) {
    fcst_start <- as.Date(fcst_start)
  }
  if (object$n_fcst > 0) {
    if (!inherits(fcst_start, "Date")) {
      tmp <- tryCatch(lubridate::ymd(rownames(object$Y)[nrow(object$Y)]), warning = function(cond) cond)
      if (inherits(tmp, "warning")) {
        stop("To summarize the forecasts, either fcst_start must be supplied or the rownames of Y be dates (YYYY-MM-DD).")
      } else {
        fcst_start <- lubridate::ymd(rownames(object$Y)[nrow(object$Y)]) %m+% months(1)
      }
    }
  }

  final_non_na <- min(unlist(apply(object$Y, 2, function(x) Position(is.na, x, nomatch = nrow(object$Y))))[object$mfbvar_prior$freq == "m"])
  final_fcst <- object$n_lags- (nrow(object$Y)-final_non_na)+1
  if (final_fcst >= 1) {
    incl_fcst <- final_fcst:(object$n_lags + object$n_fcst)
  } else {
    incl_fcst <- 1:(object$n_lags + object$n_fcst)
  }


  ret_names <- fcst_start %m+% months((-(length(incl_fcst)-object$n_fcst)):(object$n_fcst-1))
  fcst_collapsed <- tibble(variable = rep(rep(object$names_col, each = length(incl_fcst)), object$n_reps),
                           iter = rep(1:object$n_reps, each = object$n_vars*length(incl_fcst)),
                           fcst = c(object$Z_fcst[incl_fcst,,]),
                           fcst_date = rep(as.Date(as.character(ret_names)), object$n_vars*object$n_reps),
                           freq = rep(rep(object$mfbvar_prior$freq, each = length(incl_fcst)), object$n_reps),
                           time = rep(nrow(object$Y)+object$n_fcst-max(incl_fcst)+incl_fcst, object$n_vars*object$n_reps)
                           ) %>%
    transmute(variable = variable,
              iter = iter,
              year = year(fcst_date),
              quarter = quarter(fcst_date),
              fcst_date = fcst_date,
              fcst = fcst,
              freq = freq,
              time = time)
  if (aggregate_fcst) {
    fcst_collapsed <- dplyr::filter(fcst_collapsed, freq == "q") %>%
      group_by(variable, iter, year, quarter) %>%
      mutate(quarter_size = n()) %>%
      dplyr::filter(quarter_size == 3) %>%
      summarize(fcst_date = max(fcst_date), fcst = mean(fcst), freq = unique(freq), time = max(time)) %>%
      bind_rows(dplyr::filter(fcst_collapsed, freq == "m")) %>%
      ungroup()
  }

  if (!is.null(pred_bands) && !is.na(pred_bands)) {
    pred_quantiles <- c(0.5-pred_bands/2, 0.5, 0.5+pred_bands/2)
    fcst_collapsed <- group_by(fcst_collapsed, variable, time, fcst_date) %>%
      summarize(lower = quantile(fcst, prob = pred_quantiles[1]),
                median = quantile(fcst, prob = pred_quantiles[2]),
                upper = quantile(fcst, prob = pred_quantiles[3])) %>%
      ungroup()
  }


  return(fcst_collapsed)
}

#' Plot method for class \code{mfbvar_prior}
#'
#' Method for plotting \code{mfbvar_prior} objects.
#'
#' @param x object of class \code{mfbvar_prior}
#' @param nrow_facet number of rows in facet
#' @param ... Currently not in use.
#' @details The function plots the data. If the prior moments for the steady-state parameters are available in \code{x}, these are included.
#' @examples
#' prior_obj <- set_prior(Y = mf_sweden[, 4:5], freq = c("m", "q"),
#'                        n_lags = 4, n_burnin = 20, n_reps = 20, n_fcst = 4)
#' plot(prior_obj)
plot.mfbvar_prior <- function(x, nrow_facet = NULL, ...){


  ss_level <- c(0.025, 0.975)

  names_col <- if (is.null(colnames(x$Y))) paste0("x", 1:x$n_vars) else colnames(x$Y)

  if (!is.null(x$d)& !is.null(x$prior_psi_mean) & !is.null(x$prior_psi_Omega)) {
    ss_flag <- TRUE
  } else {
    ss_flag <- FALSE
  }

  if (ss_flag) {
    n_determ <- ncol(x$d)
    ss_lower  <- x$d %*% t(matrix(qnorm(ss_level[1], x$prior_psi_mean, diag(x$prior_psi_Omega)), ncol = n_determ))
    ss_median <- x$d %*% t(matrix(qnorm(0.5, x$prior_psi_mean, diag(x$prior_psi_Omega)), ncol = n_determ))
    ss_upper  <- x$d %*% t(matrix(qnorm(ss_level[2], x$prior_psi_mean, diag(x$prior_psi_Omega)), ncol = n_determ))

    ss <- data.frame(expand.grid(time = as.Date(rownames(x$Y)), variable = names_col), lower = c(ss_lower), median = c(ss_median),
                     upper = c(ss_upper))
  }

  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    row_names <- 1:nrow(x$Y)
  }
  plot_df <- data.frame(expand.grid(time = row_names, variable = names_col))
  plot_df$value <- c(as.matrix(x$Y))
  plot_df <- na.omit(plot_df)

  p <-  ggplot(mapping = aes(x = time))

  if (ss_flag) {
    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "#bdbdbd"),  alpha = 1) +
      scale_fill_manual(values = c("#bdbdbd"),
                        label = c(paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
      labs(fill = "Intervals",
           title = "Prior steady-state intervals")+
      guides(fill = guide_legend(override.aes = list(fill = c("#bdbdbd"),
                                                     linetype = c("blank"))))
  }

  p <- p +
    geom_line(data=plot_df,aes(y=value), alpha = 0.75) +
    labs(y = "Value",
         x = "Time")

  if (!ss_flag) {
    p <- p +
      labs(title = "Data")
  }

  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position="bottom")

  return(p)
}

