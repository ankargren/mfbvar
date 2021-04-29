#' Set priors for mfbvar
#'
#' The function creates an object storing all information needed for estimating a mixed-frequency BVAR. The object includes data as well as details for the model and its priors.
#'
#' @param Y data input. For monthly-quarterly data, should be a list with components containing regularly spaced time series (that inherit from \code{ts} or \code{zooreg}). If a component contains a single time series, the component itself must be named. If a component contains multiple time series, each time series must be named. Monthly variables can only contain missing values at the end of the sample, and should precede quarterly variables in the list. Matrices in which quarterly variables are padded with \code{NA} and observations stored at the end of each quarter are also accepted, but then the frequency of each variable must be given in the argument \code{freq}. Weekly-monthly mixes can be provided using the matrix way, see examples.
#' @param aggregation the aggregation scheme used for relating latent high-frequency series to their low-frequency observations. The default is \code{"average"} for averaging within each low-frequency period (e.g., quarterly observations are averages of the constituent monthly observations). The alternative \code{"triangular"} can be used for monthly-quarterly mixes, and uses the Mariano-Murasawa triangular set of weights. See details for more information.
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @param lambda2 (Only if \code{variance} is one of \code{c("diffuse", "fsv")} The cross-variable tightness
#' @templateVar lambda3 TRUE
#' @param lambda4 (Minnesota only) Prior variance of the intercept.
#' @param block_exo (Only if \code{variance} is one of \code{c("diffuse", "fsv")}) Vector of indexes/names of variables to be treated as block exogenous
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @param n_thin Store every \code{n_thin}th draw
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @param d (Steady state only) Either a matrix with same number of rows as \code{Y} and \code{n_determ} number of columns containing the deterministic terms or a string \code{"intercept"} for requesting an intercept as the only deterministic
#' term.
#' @templateVar freq TRUE
#' @param d_fcst (Steady state only) The deterministic terms for the forecasting period (not used if \code{d = "intercept"}).
#' @param prior_psi_mean (Steady state only) Vector of length \code{n_determ*n_vars} with the prior means of the steady-state parameters.
#' @param prior_psi_Omega (Steady state only) Matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)} with the prior covariance of the steady-state parameters.#'
#' @templateVar check_roots TRUE
#' @param s (Hierarchical steady state only) scalar giving the tuning parameter for the Metropolis-Hastings proposal for the kurtosis parameter. If \code{s < 0}, then adaptive Metropolis-Hastings targeting an acceptance rate of 0.44 is used, where the scaling factor is restricted to the interval \code{[-abs(s), abs(s)]}
#' @param prior_ng (Hierarchical steady state only) vector with two elements giving the parameters \code{c(c0, c1)} of the hyperprior for the global shrinkage parameter
#' @param prior_phi (Only used with common stochastic volatility) Vector with two elements \code{c(mean, variance)} for the AR(1) parameter in the log-volatility regression
#' @param prior_sigma2 (Only used with common stochastic volatility) Vector with two elements \code{c(mean, df)} for the innovation variance of the log-volatility regression
#' @param n_fac (Only used with factor stochastic volatility) Number of factors to use for the factor stochastic volatility model
#' @param n_cores (Only used with factor stochastic volatility) Number of cores to use for drawing regression parameters in parallel
#' @param ... (Only used with factor stochastic volatility) Arguments to pass along to \code{\link[factorstochvol]{fsvsample}}. See details.
#' @templateVar verbose TRUE
#' @template man_template
#' @details Some support is provided for single-frequency data sets, where \code{Y} contains variables sampled with the same frequency.
#'
#' The aggregation weights that can be used for \code{aggregation} are intra-quarterly averages (\code{aggregation = "average"}), where the quarterly observations \eqn{y_{q,t}} are assumed to relate to the underlying monthly series \eqn{z_{q,,t}} through:
#' \deqn{y_{q,t} = \frac{1}{3}(z_{q,,t} + z_{q,,t-1} + z_{q,, t-2})}
#'
#' If \code{aggregation = "triangular"}, then instead
#' \deqn{y_{q,t} = \frac{1}{9}(z_{q,,t} + 2z_{q,,t-1} + 3z_{q,, t-2}) + 2z_{q,, t-3}) + z_{q,, t-4})}
#'
#' The latter is typically used when modeling growth rates, and the former when working with log-levels.
#'
#' If the steady-state prior is to be used, the deterministic matrix needs to be supplied, or a string indicating that the intercept should be the only deterministic term (\code{d = "intercept"}). If the latter, \code{d_fcst} is automatically set to be intercept only. Otherwise, if forecasts are requested
#' (\code{n_fcst > 0}) also \code{d_fcst} must be provided. Finally, the prior means of the steady-state parameters must (at the very minimum) also be
#' provided in \code{prior_psi_mean}. The steady-state prior involves inverting the lag polynomial. For this reason, draws in which the largest eigenvalue
#' (in absolute value) of the lag polynomial is greater than 1 are discarded and new draws are made if \code{check_roots = TRUE}. The maximum number of
#' attempts is 1,000.
#'
#' For modeling stochastic volatility by the factor stochastic volatility model, the number of factors to use must be supplied. Further arguments can be passed along, but are not included as formal arguments. If the default settings are not overriden, the defaults used are as follows (see \code{\link[factorstochvol]{fsvsample}} for descriptions):
#' \itemize{
#'   \item{\code{priormu}}{\code{ = c(0, 10)}}
#'   \item{\code{priorphiidi}}{\code{ = c(10, 3)}}
#'   \item{\code{priorphifac}}{\code{ = c(10, 3)}}
#'   \item{\code{priorsigmaidi}}{\code{ = 1}}
#'   \item{\code{priorsigmafac}}{\code{ = 1}}
#'   \item{\code{priorfacload}}{\code{ = 1}}
#'   \item{\code{restrict}}{\code{ = "none"}}
#' }
#'
#'
#' @return An object of class \code{mfbvar_prior} that is used as input to \code{estimate_mfbvar}.
#' @examples
#' # Standard list-based way
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' prior_obj <- set_prior(prior_obj, n_fcst = 4)
#'
#' # Weekly-monthly mix of data, four weeks per month
#' Y <- matrix(rnorm(400), 100, 4)
#' Y[setdiff(1:100,seq(4, 100, by = 4)), 4] <- NA
#' prior_obj <- set_prior(Y = Y, freq = c(rep("w", 3), "m"),
#'                        n_lags = 4, n_reps = 10)
#' @seealso \code{\link{estimate_mfbvar}}, \code{\link{interval_to_moments}}, \code{\link{print.mfbvar_prior}}, \code{\link{summary.mfbvar_prior}}, \code{\link[factorstochvol]{fsvsample}}
set_prior <- function(prior_obj, Y, aggregation, prior_Pi_AR1, lambda1,
                      lambda2, lambda3, lambda4, block_exo, n_lags, n_fcst,
                      n_thin, n_reps, n_burnin, freq, d, d_fcst, prior_psi_mean,
                      prior_psi_Omega, check_roots, s, prior_ng, prior_phi,
                      prior_sigma2, n_fac, n_cores, priormu, priorphiidi,
                      priorphifac, priorsigmaidi, priorsigmafac, priorfacload,
                      restrict, verbose, ...) {
  if (!missing(prior_obj)) {
    prior_call <- mget(ls())
    prior_call <- prior_call[!(names(prior_call) != "prior_obj")]
    call_names <- names(prior_call)
    prior_obj <- prior_obj[which(!(names(prior_obj) %in% call_names))]
    prior_obj <- c(prior_obj, prior_call)
    prior_obj$supplied_args <- union(prior_obj$supplied_args, call_names)
    ret <- check_prior(prior_obj)
  } else {
    prior_call <- mget(ls())
    prior_call <- prior_call[which(names(prior_call) %in% names(as.list(match.call()))[-1])]
    prior_call$supplied_args <- names(prior_call)
    ret <- check_prior(prior_call)
  }
  class(ret) <- "mfbvar_prior"
  ret$call <- match.call()
  return(ret)
}

call_wrapper <- function(prior_call, defaults, prev_args) {
  if ("prior_obj" %in% names(prior_call)) {
    prior_obj <- prior_call[["prior_obj"]]
    prev_args <- prior_obj$wrapper_args
    prior_call <- prior_call[!(names(prior_call) == "prior_obj")]
    prior_obj <- prior_obj[!(names(prior_obj) %in% names(prior_call))]
    prior_call <- c(prior_obj, prior_call)
  } else {
    prev_args <- NULL
  }
  supplied_args <- names(prior_call)
  defaults <- defaults[!(names(defaults) %in% supplied_args)]
  prior_call <- append(prior_call,
                       defaults)
  prior_obj <- do.call("set_prior", prior_call)
  prior_obj$wrapper_args <- supplied_args
  return(prior_obj)
}

#' @rdname set_prior
set_init <- function(prior_obj, Y, aggregation = "average", n_lags, n_fcst = 0,
                     n_reps, n_burnin, n_thin = 1, freq = NULL,
                     verbose = FALSE) {
  prior_call <- mget(ls())
  prior_call <- prior_call[which(names(prior_call) %in% names(as.list(match.call()))[-1])]
  defaults <- formals()
  prev_args <- ifelse(missing(prior_obj), "", prior_obj$wrapper_args)
  if (missing(prior_obj)) {
    prev_args <- NULL
  } else {
    prev_args <- prior_obj$wrapper_args
  }
  prior_obj <- call_wrapper(prior_call, defaults, prev_args)
  return(prior_obj)
}

#' @rdname set_prior
set_prior_minn <- function(prior_obj, prior_Pi_AR1 = 0, lambda1 = 0.2,
                           lambda2 = 0.5, lambda3 = 1, lambda4 = 10000,
                           block_exo = NULL, check_roots = FALSE) {
  prior_call <- mget(ls())
  defaults <- formals()
  if (!missing(prior_obj)) {
    prior_call <- prior_call[which(names(prior_call) %in% names(as.list(match.call()))[-1])]
  }
  prior_obj <- call_wrapper(prior_call, defaults)
  return(prior_obj)
}

#' @rdname set_prior
set_prior_ss <- function(prior_obj, prior_psi_mean, prior_psi_Omega,
                         d = "intercept", d_fcst = NULL) {
  prior_call <- mget(ls())
  defaults <- formals()
  prev_args <- ifelse(missing(prior_obj), "", prior_obj$wrapper_args)
  if (missing(prior_obj)) {
    prior_call <- prior_call[-1]
    prev_args <- NULL
  } else {
    prev_args <- prior_obj$wrapper_args
  }
  prior_obj <- call_wrapper(prior_call, defaults, prev_args)
  return(prior_obj)
}

#' @rdname set_prior
set_prior_ssng <- function(prior_obj, s = -1000, prior_ng = c(0.01, 0.01)) {
  prior_call <- mget(ls())
  defaults <- formals()
  prev_args <- ifelse(missing(prior_obj), "", prior_obj$wrapper_args)
  if (missing(prior_obj)) {
    prior_call <- prior_call[-1]
    prev_args <- NULL
  } else {
    prev_args <- prior_obj$wrapper_args
  }
  prior_obj <- call_wrapper(prior_call, defaults, prev_args)
  return(prior_obj)
}

#' @rdname set_prior
set_prior_fsv <- function(prior_obj, n_fac, n_cores = 1, priormu = c(0, 10),
                          priorphiidi = c(10, 3), priorphifac = c(10, 3),
                          priorsigmaidi = 1, priorsigmafac = 1,
                          priorfacload = 1, restrict = "none") {
  prior_call <- mget(ls())
  defaults <- formals()
  prev_args <- ifelse(missing(prior_obj), "", prior_obj$wrapper_args)
  if (missing(prior_obj)) {
    prior_call <- prior_call[-1]
    prev_args <- NULL
  } else {
    prev_args <- prior_obj$wrapper_args
  }
  prior_obj <- call_wrapper(prior_call, defaults, prev_args)
  return(prior_obj)
}

#' @rdname set_prior
set_prior_csv <- function(prior_obj, prior_phi = c(0.9, 0.1),
                          prior_sigma2 = c(0.01, 4)) {
  prior_call <- mget(ls())
  defaults <- formals()
  prev_args <- ifelse(missing(prior_obj), "", prior_obj$wrapper_args)
  if (missing(prior_obj)) {
    prior_call <- prior_call[-1]
    prev_args <- NULL
  } else {
    prev_args <- prior_obj$wrapper_args
  }
  prior_obj <- call_wrapper(prior_call, defaults, prev_args)
  return(prior_obj)
}


check_prior <- function(prior_obj) {
  is_numeric_vector <- function(x, n, null_ok = TRUE, matrix_ok = FALSE) {
    null <- (null_ok && is.null(x))
    num <- is.numeric(x)
    vec <- is.vector(x)
    mat <- (is.matrix(x) && min(dim(x)) == 1 && matrix_ok)
    if (missing(n)) {
      null || (num && (vec || mat))
    } else {
      null || (num && (vec || mat) && length(x) %in% n)
    }
  }

  is_scalar <- function(x, null_ok = TRUE) {
    is_numeric_vector(x, n = 1, null_ok = null_ok, matrix_ok = FALSE)
  }

  assert_class <- function(prior_obj, x, y) {
    if (!inherits(prior_obj[[x]], y)) {
      stop(sprintf("%s must inherit from class %s.",
                   x, paste(y, collapse = " or ")))
    }
  }

  assert_numeric_vector <- function(prior_obj, x, n = 2, null_ok = TRUE,
                                    matrix_ok = FALSE) {
    if (!is_numeric_vector(prior_obj[[x]], n, null_ok, matrix_ok)) {
      stop(sprintf("%s must be a %d-dimensional numeric vector.", x, n))
    }
  }

  assert_scalar <- function(prior_obj, x, null_ok = FALSE) {
    assert_numeric_vector(prior_obj, x, n = 1, null_ok = null_ok,
                          matrix_ok = FALSE)
  }

  scalar_values <- c("lambda1", "lambda2", "lambda3", "lambda4", "s", "n_lags",
                     "n_thin", "n_reps", "n_burnin", "a", "n_fcst")

  required <- c("Y", "n_lags")
  missing_required <- !(required %in% names(prior_obj))
  if (any(missing_required)) {
    stop(sprintf("Missing required arguments: %s",
                 paste(required[missing_required],
                       collapse = ", ")))
  }

  assert_class(prior_obj, "Y", c("matrix", "data.frame", "list"))
  if (!is.matrix(prior_obj$Y)) {
    if (inherits(prior_obj$Y, "list")) {
      list_conv <- list_to_matrix(prior_obj$Y)
      prior_obj$Y <- list_conv[[1]]
      prior_obj$freq <- list_conv[[2]]
    } else if (is.data.frame(prior_obj$Y)) {
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

  if (nrow(prior_obj$Y) >
      max(unlist(apply(prior_obj$Y, 2, function(x)
        Position(is.na, x, nomatch = nrow(prior_obj$Y)))))) {
    stop("Y: remove final rows containing only NAs.")
  }


  if ("freq" %in% prior_obj$supplied_args) {
    assert_class(prior_obj, "freq", "character")
    if (!all(prior_obj$freq %in% c("w", "m", "q"))) {
      stop("Elements of freq must be 'w', 'm' or 'q'.")
    } else if (length(prior_obj$freq) != ncol(prior_obj$Y)) {
      stop("The length of freq is ", length(prior_obj$freq), ", but Y has ", ncol(prior_obj$Y), " columns.")
    } else {
      freq_pos <- c(
        ifelse(any(prior_obj$freq == "q"), which.max(prior_obj$freq == "q"), NA),
        ifelse(any(prior_obj$freq == "m"), which.max(prior_obj$freq == "m"), NA),
        ifelse(any(prior_obj$freq == "w"), which.max(prior_obj$freq == "w"), NA)
      )
      freqs <- c("q", "m", "w")
      freqs <- freqs[!is.na(freq_pos)]
      if (length(freqs) > 2) {
        stop("mfbvar can currently only handle a mix of two frequencies.")
      }
      if (length(freqs)>1 && freqs[1]=="q" && freqs[2] == "w") {
        stop("mfbvar can currently only handle weekly-monthly or monthly-quarterly mixes.")
      }
      prior_obj$freqs <- freqs
      if (length(freq_pos[!is.na(freq_pos)])>1 && diff(freq_pos[!is.na(freq_pos)])>0) {
        stop("Variables must be placed in weekly-monthly-quarterly order.")
      }
    }
  }

  if (length(unique(prior_obj$freq)) > 1) {
    if (min(unlist(
      apply(prior_obj$Y[, prior_obj$freq %in% freqs[-1], drop = FALSE], 2,
            function(x) Position(is.na, x, nomatch = 1e10)))) == 1) {
      stop("Y: high-frequency variables are NA at the beginning of the sample.")
    }
  } else {
    print(str(prior_obj$Y))
    if (any(is.na(prior_obj$Y))) {
      stop("Y: single-frequency estimation requires the data to contain no NAs.")
    }
  }

  intercept_flag <- FALSE

  if ("d" %in% prior_obj$supplied_args) {
    if (is_scalar(prior_obj$d) || prior_obj$d == "intercept") {
      intercept_flag <- TRUE
      prior_obj$d <- matrix(1, nrow = nrow(prior_obj$Y), 1)
    } else if (all(prior_obj$d == 1)) {
      intercept_flag <- TRUE
    }

    prior_obj$intercept_flag <- intercept_flag

    if (!intercept_flag) {
      if ("d_fcst" %in% prior_obj$supplied_args) {
        assert_class(prior_obj, "d_fcst", "matrix")
        if (ncol(prior_obj$d) != ncol(prior_obj$d_fcst)) {
          stop("d has", ncol(prior_obj$d), " columns and d_fcst", ncol(prior_obj$d_fcst), ".")
        }
      }
    }

    if (nrow(prior_obj$Y) != nrow(prior_obj$d)) {
      stop("Y has ", nrow(prior_obj$Y), "rows and d ", nrow(prior_obj$d), "rows, but they must be equal.")
    }

    assert_class(prior_obj, "d", "matrix")
  }

  assert_class(prior_obj, "aggregation", c("character", "matrix"))
  if (is.matrix(prior_obj$aggregation)) {
    prior_obj$Lambda <- prior_obj$aggregation
    prior_obj$aggregation <- "custom"
  } else {
    freq <- prior_obj$freq
    freqs <- prior_obj$freqs
    n_l <- ifelse(length(freqs)>1, sum(freq == freqs[1]), 0)
    n_h <- ifelse(length(freqs)>1, sum(freq == freqs[2]), length(freq))
    if (n_l > 0 && freqs[1] == "q") {
      if (prior_obj$aggregation == "average") {
        prior_obj$Lambda_ <- build_Lambda(rep("average", n_l), 3)
      } else {
        prior_obj$Lambda_ <- build_Lambda(rep("triangular", n_l), 5)}
    } else if (n_l == 0) {
      prior_obj$Lambda_ <- diag(ncol(prior_obj$Y))
    } else if (freqs[1] == "m") {
      if (prior_obj$aggregation == "triangular") {
        stop("Triangular aggregation not supported for weekly data.")
      } else {
        prior_obj$Lambda_ <- matrix(0.25, 1, 4)
        prior_obj$Lambda_ <- kronecker(prior_obj$Lambda_, diag(n_l))
      }
    }
  }

  if (prior_obj$aggregation == "triangular") {
    if (prior_obj$n_lags < 5) {
      stop("The number of lags must be at least 5 when using triangular aggregation.")
    }
  } else if (prior_obj$aggregation == "average") {
    if (prior_obj$n_lags < 3) {
      stop("The number of lags must be at least 3 when using intra-quarterly averaging.")
    }
  }

  assert_numeric_vector(prior_obj, "prior_Pi_AR1", c(1, ncol(prior_obj$Y)))
  if (is_scalar(prior_obj$prior_Pi_AR1)) {
    prior_obj$prior_Pi_AR1 <- rep(prior_obj$prior_Pi_AR1, ncol(prior_obj$Y))
  }

  for (i in seq_along(scalar_values)) {
    if (scalar_values[i] %in% prior_obj$supplied_args) {
      assert_scalar(prior_obj, scalar_values[i], null_ok = TRUE)
    }
  }

  if ("block_exo" %in% prior_obj$supplied_args) {
    if (!is.null(prior_obj$block_exo) && !is.vector(prior_obj$block_exo)) {
      stop("block_exo must be a vector of indexes or names.")
    } else {
      if (is.character(prior_obj$block_exo)) {
        if (all(prior_obj$block_exo %in% colnames(prior_obj$Y))) {
          prior_obj$block_exo <- which(prior_obj$block_exo %in% colnames(prior_obj$Y))
        }
      }
    }
  }

  if (!is_numeric_vector(prior_obj$prior_psi_mean)) {
    stop("prior_psi_mean must be a vector or matrix with one row or column.")
  } else {
    prior_obj$prior_psi_mean <- c(prior_obj$prior_psi_mean)
    if (length(prior_obj$prior_psi_mean) %% ncol(prior_obj$Y) != 0) {
      stop("prior_psi_mean has ", length(prior_obj$prior_psi_mean), " elements,
           but there are ", ncol(prior_obj$Y), " variables in Y.")
    }
  }

  if ("prior_psi_Omega" %in% prior_obj$supplied_args) {
    assert_class(prior_obj, "prior_psi_Omega", "matrix")
    if (dim(prior_obj$prior_psi_Omega)[1] != dim(prior_obj$prior_psi_Omega)[2]) {
      stop("prior_psi_Omega must be a positive-definite symmetric matrix.")
    }
    if (dim(prior_obj$prior_psi_Omega)[1] != length(prior_obj$prior_psi_mean)) {
      stop("The dimension of prior_psi_Omega must correspond to the number of elements in prior_psi_mean.")
    }
  }

  if (intercept_flag && prior_obj$n_fcst > 0) {
    prior_obj$d_fcst <- matrix(1, nrow = prior_obj$n_fcst, ncol = 1)
  }

  if (!is.null(prior_obj$check_roots) && !is.logical(prior_obj$check_roots)) {
    stop("check_roots: must be logical.")
  }

  assert_numeric_vector(prior_obj, "prior_ng", c(1, 2), null_ok = TRUE)
  assert_numeric_vector(prior_obj, "prior_phi", null_ok = TRUE)
  assert_numeric_vector(prior_obj, "prior_sigma2", null_ok = TRUE)
  assert_scalar(prior_obj, "n_fac", null_ok = TRUE)
  assert_scalar(prior_obj, "n_cores", null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priormu", n = 2, null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priorphiidi", n = 2, null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priorphifac", n = 2, null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priorsigmaidi", n = c(1, ncol(prior_obj$Y)),
                        null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priorsigmafac", n = c(1, prior_obj$n_fac),
                        null_ok = TRUE)
  assert_numeric_vector(prior_obj, "priorsigmafac", n = c(1, prior_obj$n_fac),
                        null_ok = TRUE)
  if (!is.null(prior_obj$priorfacload)) {
    if (!(is.numeric(prior_obj$priorfacload) && (length(prior_obj$priorfacload) == 1 || dim(prior_obj$supplied_args$priorfacload) == c(ncol(prior_obj$Y), prior_obj$factors)))) {
      stop(sprintf("priorfacload should be a scalar value or an n_vars x n_fac matrix, but is %s with %d elements", class(prior_obj$priorfacload), length(prior_obj$priorfacload)))
    }
  }

  if (!is.null(prior_obj$restrict)) {
    if (!(is.character(prior_obj$restrict) && length(prior_obj$restrict) == 1)) {
      stop(sprintf("restrict should be a single string, but is %s of length %d", class(prior_obj$restrict), length(prior_obj$restrict)))
    } else {
      if (!(prior_obj$restrict %in% c("none", "upper"))) {
        stop(sprintf("restrict should be 'none' or 'upper', but is %s", prior_obj$restrict))
      }
    }
  }
  return(prior_obj)
}


#' Mixed-frequency Bayesian VAR
#'
#' The main function for estimating a mixed-frequency BVAR.
#'
#' @param mfbvar_prior a \code{mfbvar_prior} object
#' @param prior either \code{"ss"} (steady-state prior), \code{"ssng"} (hierarchical steady-state prior with normal-gamma shrinkage) or \code{"minn"} (Minnesota prior)
#' @param variance form of the error variance-covariance matrix: \code{"iw"} for the inverse Wishart prior, \code{"diffuse"} for a diffuse prior, \code{"csv"} for common stochastic volatility or \code{"fsv"} for factor stochastic volatility
#' @param ... additional arguments to \code{update_prior} (if \code{mfbvar_prior} is \code{NULL}, the arguments are passed on to \code{set_prior})
#' @return
#'  An object of class \code{mfbvar}, \code{mfbvar_<prior>} and \code{mfbvar_<prior>_<variance>} containing posterior quantities as well as the prior object. For all choices of \code{prior} and \code{variance}, the returned object contains:
#' \item{Pi}{Array of dynamic coefficient matrices; \code{Pi[,, r]} is the \code{r}th draw}
#' \item{Z}{Array of monthly processes; \code{Z[,, r]} is the \code{r}th draw}
#' \item{Z_fcst}{Array of monthly forecasts; \code{Z_fcst[,, r]} is the \code{r}th forecast. The first \code{n_lags}
#' rows are taken from the data to offer a bridge between observations and forecasts and for computing nowcasts (i.e. with ragged edges).}
#' \subsection{Steady-state priors}{
#' If \code{prior = "ss"}, it also includes:
#' \describe{\item{\code{psi}}{Matrix of steady-state parameter vectors; \code{psi[r,]} is the \code{r}th draw}
#' \item{\code{roots}}{The maximum eigenvalue of the lag polynomial (if \code{check_roots = TRUE})}}
#'
#' If \code{prior = "ssng"}, it also includes:
#' \describe{
#' \item{\code{psi}}{Matrix of steady-state parameter vectors; \code{psi[r,]} is the \code{r}th draw}
#' \item{\code{roots}}{The maximum eigenvalue of the lag polynomial (if \code{check_roots = TRUE})}
#' \item{\code{lambda_psi}}{Vector of draws of the global hyperparameter in the normal-Gamma prior}
#' \item{\code{phi_psi}}{Vector of draws of the auxiliary hyperparameter in the normal-Gamma prior}
#' \item{\code{omega_psi}}{Matrix of draws of the prior variances of psi; \code{omega_psi[r, ]} is the \code{r}th draw, where \code{diag(omega_psi[r, ])} is used as the prior covariance matrix for psi}}}
#' \subsection{Constant error covariances}{
#' If \code{variance = "iw"} or \code{variance = "diffuse"}, it also includes:
#' \describe{\item{\code{Sigma}}{Array of error covariance matrices; \code{Sigma[,, r]} is the \code{r}th draw}}}
#' \subsection{Time-varying error covariances}{
#' If \code{variance = "csv"}, it also includes:
#' \describe{\item{\code{Sigma}}{Array of error covariance matrices; \code{Sigma[,, r]} is the \code{r}th draw}
#' \item{\code{phi}}{Vector of AR(1) parameters for the log-volatility regression; \code{phi[r]} is the \code{r}th draw}
#' \item{\code{sigma}}{Vector of error standard deviations for the log-volatility regression; \code{sigma[r]} is the \code{r}th draw}
#' \item{\code{f}}{Matrix of log-volatilities; \code{f[r, ]} is the \code{r}th draw}}
#'
#' If \code{variance = "fsv"}, it also includes:
#' \describe{\item{\code{facload}}{Array of factor loadings; \code{facload[,, r]} is the \code{r}th draw}
#' \item{\code{latent}}{Array of latent log-volatilities; \code{latent[,, r]} is the \code{r}th draw}
#' \item{\code{mu}}{Matrix of means of the log-volatilities; \code{mu[, r]} is the \code{r}th draw}
#' \item{\code{phi}}{Matrix of AR(1) parameters for the log-volatilities; \code{phi[, r]} is the \code{r}th draw}
#' \item{\code{sigma}}{Matrix of innovation variances for the log-volatilities; \code{sigma[, r]} is the \code{r}th draw}}}
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{predict.mfbvar}}, \code{\link{plot.mfbvar_minn}},
#' \code{\link{plot.mfbvar_ss}}, \code{\link{varplot}}, \code{\link{summary.mfbvar}}
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' @references
#' Ankargren, S., Unosson, M., & Yang, Y. (2020) A Flexible Mixed-Frequency Bayesian Vector Autoregression with a Steady-State Prior. \emph{Journal of Time Series Econometrics}, 12(2), \doi{10.1515/jtse-2018-0034}.\cr
#' Ankargren, S., & Jonéus, P. (2020) Simulation Smoothing for Nowcasting with Large Mixed-Frequency VARs. \emph{Econometrics and Statistics}, \doi{10.1016/j.ecosta.2020.05.007}.\cr
#' Ankargren, S., & Jonéus, P. (2019) Estimating Large Mixed-Frequency Bayesian VAR Models. arXiv:1912.02231, \url{https://arxiv.org/abs/1912.02231}.\cr
#' Kastner, G., & Huber, F. (2020) Sparse Bayesian Vector Autoregressions in Huge Dimensions. \emph{Journal of Forecasting}, 39, 1142--1165. \doi{10.1002/for.2680}.\cr
#' Schorfheide, F., & Song, D. (2015) Real-Time Forecasting With a Mixed-Frequency VAR. \emph{Journal of Business & Economic Statistics}, 33(3), 366--380. \doi{10.1080/07350015.2014.954707}\cr

estimate_mfbvar <- function(mfbvar_prior = NULL, prior, variance = "iw", ...) {
  time_out <- Sys.time()
  args <- list(...)
  if (hasArg(mfbvar_prior)) {
    if (!inherits(mfbvar_prior, "mfbvar_prior")) {
      stop("mfbvar_prior must be of class mfbvar_prior.")
    } else {
      if (length(args) > 0) {
        mfbvar_prior <- set_prior(mfbvar_prior, ...)
      }
    }
  } else {
    mfbvar_prior <- set_prior(...)
  }

  if (hasArg(prior_type)) {
    warning("The argument 'prior_type' is deprecated (starting in 0.5.0). Use 'prior' instead.", call. = FALSE)
    prior <- args$prior_type
  }

  prior_opts <- c("ss", "ssng", "minn", "dl")
  variance_opts <- c("iw", "fsv", "csv", "diffuse")
  envir <- environment()
  if (!(prior %in% prior_opts)) {
    stop(sprintf("prior must be %s or %s.",
                 paste(prior_opts[-length(prior_opts)], collapse = ", "),
                 prior_opts[length(prior_opts)]))
  } else {
    assign(prior, TRUE, envir)
    vapply(setdiff(prior_opts, prior),
           function(x) assign(x, FALSE, envir),
           logical(1))
  }
  if (!(variance %in% variance_opts)) {
    stop(sprintf("variance must be %s or %s.",
                 paste(variance_opts[-length(variance_opts)], collapse = ", "),
                 variance_opts[length(variance_opts)]))
  } else {
    assign(variance, TRUE, envir)
    vapply(setdiff(variance_opts, variance),
           function(x) assign(x, FALSE, envir),
           logical(1))
  }

  if (prior == "dl" && !(variance %in% c("fsv", "diffuse"))) {
    stop("The Dirichlet-Laplace prior (dl) can only be used with
         variance specifications fsv and diffuse.")
  }

  class(mfbvar_prior) <- c(sprintf("mfbvar_%s_%s", prior, variance),
                           sprintf("mfbvar_%s", prior),
                           sprintf("mfbvar_%s", variance),
                           class(mfbvar_prior))

  time_out <- c(time_out, Sys.time())
  main_run <- mfbvar_sampler(mfbvar_prior,
                             minn = minn,
                             ss = ss,
                             ssng = ssng,
                             dl = dl,
                             iw = iw,
                             csv = csv,
                             diffuse = diffuse,
                             fsv = fsv,
                             ...)
  time_out <- c(time_out, Sys.time())
  if (mfbvar_prior$verbose) {
    time_diff <- Sys.time() - time_out[1]
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
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(length(names_row)-main_run$n_lags+1):length(names_row)]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col
  } else {
    names_fcst <- NULL
  }


  main_run$names_row <- names_row
  main_run$names_col <- names_col
  main_run$names_fcst <- names_fcst
  main_run$mfbvar_prior <- mfbvar_prior

  dimnames(main_run$Z) <- list(time = names_row[(nrow(mfbvar_prior$Y)-nrow(main_run$Z)+1):nrow(mfbvar_prior$Y)],
                               variable = names_col,
                               iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))

  if (variance %in% c("iw", "diffuse")) {
    dimnames(main_run$Sigma) <- list(names_col,
                                     names_col,
                                     iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
  }


  if (prior %in% c("ss", "ssng")) {
    if (is.null(colnames(mfbvar_prior$d))) {
      names_determ <- paste0("d", 1:ncol(mfbvar_prior$d))
    } else {
      names_determ <- colnames(mfbvar_prior$d)
    }
    rownames(mfbvar_prior$d) <- rownames(mfbvar_prior$Y)
    main_run$names_determ <- names_determ
    n_determ <- dim(mfbvar_prior$d)[2]
    dimnames(main_run$psi) <- list(param = paste0(rep(names_col, n_determ), ".", rep(names_determ, each = n_vars)),
                                   NULL,
                                   iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars)),
                                  iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
  } else {
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = c("const", paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars))),
                                  iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
  }

  if (length(mfbvar_prior$freqs) == 1) {
    class(main_run) <- c("sfbvar", sprintf("sfbvar_%s_%s", prior, variance), sprintf("sfbvar_%s", prior), sprintf("sfbvar_%s", variance))
  } else {
    class(main_run) <- c("mfbvar", sprintf("mfbvar_%s_%s", prior, variance), sprintf("mfbvar_%s", prior), sprintf("mfbvar_%s", variance))
  }
  time_out <- c(time_out, Sys.time())
  main_run$time_out <- time_out
  main_run$variance <- variance
  main_run$prior <- prior
  return(main_run)
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
  n_reps <- x$n_reps
  n_T <- x$n_T_
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
    variances_fsv(variances, x$latent, x$facload, variables_num, n_fac, n_reps, n_T, n_vars, n_plotvars)
  }
  if (sv_type == "csv") {
    variances_csv(variances, x$Sigma, x$latent, n_T, n_reps, variables_num)
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
