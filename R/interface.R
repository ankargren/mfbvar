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
#' The function \code{update_prior} can be used to update an existing prior object. See the examples.
#'
#' @return An object of class \code{mfbvar_prior} that is used as input to \code{estimate_mfbvar}.
#' @examples
#' # Standard list-based way
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' prior_obj <- update_prior(prior_obj, n_fcst = 4)
#'
#' # Weekly-monthly mix of data, four weeks per month
#' Y <- matrix(rnorm(400), 100, 4)
#' Y[setdiff(1:100,seq(4, 100, by = 4)), 4] <- NA
#' prior_obj <- set_prior(Y = Y, freq = c(rep("w", 3), "m"),
#'                        n_lags = 4, n_reps = 10)
#' @seealso \code{\link{estimate_mfbvar}}, \code{\link{update_prior}}, \code{\link{interval_to_moments}}, \code{\link{print.mfbvar_prior}}, \code{\link{summary.mfbvar_prior}}, \code{\link[factorstochvol]{fsvsample}}
set_prior <- function(Y, aggregation = "average", prior_Pi_AR1 = 0, lambda1 = 0.2,
                      lambda2 = 0.5, lambda3 = 1, lambda4 = 10000, block_exo = NULL, n_lags,
                      n_fcst = 0, n_thin = 1, n_reps, n_burnin = n_reps, freq = NULL, d = NULL, d_fcst = NULL,
                      prior_psi_mean = NULL, prior_psi_Omega = NULL, check_roots = FALSE,
                      s = -1000, prior_ng = c(0.01, 0.01),
                      prior_phi = c(0.9, 0.1),
                      prior_sigma2 = c(0.01, 4), n_fac = NULL,
                      n_cores = 1, verbose = FALSE, ...) {
  prior_call <- mget(names(formals())[names(formals()) != "..."], sys.frame(sys.nframe()))
  prior_call$supplied_args <- names(as.list(match.call()))[-1]
  ellipsis <- list(...)
  fsv_names <- names(ellipsis)
  fsv_arguments <- c("priormu", "priorphiidi", "priorphifac", "priorsigmaidi", "priorsigmafac", "priorfacload")
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
#' @param prior_obj an object of class \code{mfbvar_prior}
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

check_prior <- function(prior_obj) {
  if (!is.matrix(prior_obj$Y)) {
    if (inherits(prior_obj$Y, "list")) {
      list_conv <- list_to_matrix(prior_obj$Y)
      prior_obj$Y <- list_conv[[1]]
      prior_obj$freq <- list_conv[[2]]
      prior_obj$supplied_args <- c(prior_obj$supplied_args, "freq")
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
    } else {
      stop(paste0("Y is of class ", class(prior_obj$Y), ", but must be matrix, data frame/tibble or a list of ts or zooreg objects."))
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
    } else if (!all(prior_obj$freq %in% c("w", "m", "q"))) {
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
      if (length(freqs)>2) {
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
  } else {
    stop("freq: must be supplied.")
  }

  if (length(freqs)>1) {
      if (min(unlist(apply(prior_obj$Y[, prior_obj$freq %in% freqs[-1], drop = FALSE], 2, function(x) Position(is.na, x, nomatch = 9999999999)))) == 1) {
        stop("Y: high-frequency variables are NA at the beginning of the sample.")
      }
  } else {
      if (any(is.na(prior_obj$Y))) {
        stop("Y: single-frequency estimation requires the data to contain no NAs.")
      }
  }


  if ("aggregation" %in% prior_obj$supplied_args) {
    if (is.atomic(prior_obj$aggregation) || is.matrix(prior_obj$aggregation)) {
    } else {
      stop("aggregation must be a character vector or a matrix, but is now of class ", class(prior_obj$aggregation), ".")
    }
  } else {
    prior_obj$aggregation <- "average"
  }

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

  if ("block_exo" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$block_exo)) {
      stop("block_exo must be a vector of indexes or names.")
    } else {
      if (is.character(prior_obj$block_exo)) {
        if (all(prior_obj$block_exo %in% colnames(prior_obj$Y))) {
          prior_obj$block_exo <- which(prior_obj$block_exo %in% colnames(prior_obj$Y))
        }
      }
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "block_exo")
  }


  if ("prior_psi_mean" %in% prior_obj$supplied_args) {
    if (!(is.atomic(prior_obj$prior_psi_mean) || is.matrix(prior_obj$prior_psi_mean))) {
      stop("prior_psi_mean must be a vector or matrix with one row or column.")
    }
    if (is.atomic(prior_obj$prior_psi_mean)) {
      if (length(prior_obj$prior_psi_mean) %% ncol(prior_obj$Y) != 0) {
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
      if (dim(prior_obj$prior_psi_Omega)[1] != length(prior_obj$prior_psi_mean)) {
        stop("The dimension of prior_psi_Omega must correspond to the number of elements in prior_psi_mean.")
      }
    }
  }

  if ("s" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$s) || length(prior_obj$s) > 1) {
      stop("s must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "s")
  }


  if ("prior_ng" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$prior_ng) || length(prior_obj$prior_ng) > 2) {
      stop("prior_ng must be a vector with one or two elements.")
    } else {
      if (length(prior_obj$prior_ng) == 1) {
        prior_obj$prior_ng <- c(prior_obj$prior_ng, prior_obj$prior_ng)
      }
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "prior_ng")
  }

  if ("n_lags" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_lags) || length(prior_obj$n_lags) > 1) {
      stop("n_lags must be a vector with a single element.")
    }
  } else {
    stop("n_lags: No lag length specified.\n")
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

  if ("n_reps" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$n_reps) || length(prior_obj$n_reps) > 1) {
      stop("n_reps must be a vector with a single element.")
    }
  } else {
    stop("n_reps: Number of draws to use in main chain not specified.\n")
  }

  if (!is.atomic(prior_obj$n_burnin) || length(prior_obj$n_burnin) > 1) {
    stop("n_burnin must be a vector with a single element.")
  } else if (!("n_burnin" %in% prior_obj$supplied_args)) {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "n_burnin")
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

    if (!is.atomic(prior_obj$n_cores) || length(prior_obj$n_cores) > 1) {
      stop("n_cores must be a vector with a single element.")
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
      if (!(is.numeric(prior_obj$priorsigmaidi))) {
        stop("priorsigmaidi should be numeric.")
      }
      if (length(prior_obj$priorsigmaidi) == 1) {
      } else if (length(prior_obj$priorsigmaidi) == ncol(prior_obj$Y)) {
      } else {
        stop("priorsigmaidi should be a numeric vector of length 1 or n_vars.")
      }
    } else {
      prior_obj$priorsigmaidi <- rep(1, ncol(prior_obj$Y))
    }

    if ("priorsigmafac" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$priorsigmafac))) {
        stop("priorsigmafac should be numeric.")
      }
      if (length(prior_obj$priorsigmafac) == 1) {
      } else if (length(prior_obj$priorsigmafac) == prior_obj$n_fac) {
      } else {
        stop("priorsigmafac should be a numeric vector of length 1 or n_fac")
      }
    } else {
      prior_obj$priorsigmafac <- rep(1, prior_obj$n_fac)
    }

    if ("priorfacload" %in% prior_obj$supplied_args) {
      if (!(is.numeric(prior_obj$supplied_args$priorfacload) && (length(prior_obj$supplied_args$priorfacload) == 1 || dim(prior_obj$supplied_args$priorfacload) == c(ncol(prior_obj$Y), prior_obj$factors)))) {
        stop(sprintf("priorfacload should be a scalar value or an n_vars x n_fac matrix, but is %s with %d elements", class(prior_obj$priorfacload), length(prior_obj$priorfacload)))
      }
    } else {
      prior_obj$priorfacload <- 1
    }

    if ("restrict" %in% prior_obj$supplied_args) {
      if (!(is.character(prior_obj$restrict) && length(prior_obj$priorng) == 1)) {
        stop(sprintf("restrict should be a single string, but is %s of length %d", class(prior_obj$restrict), length(prior_obj$restrict)))
      } else {
        if (!(prior_obj$restrict %in% c("none", "upper"))) {
          stop(sprintf("restrict should be 'none' or 'upper', but is %s", prior_obj$restrict))
        }
      }
    } else {
      prior_obj$restrict <- "none"
    }



  } else if (is.null(prior_obj$n_fac) && any(prior_obj$supplied_args %in% c("priormu", "priorphiidi", "priorphifac", "priorsigmaidi", "priorsigmafac",
                   "priorfacload", "restrict"))) {
    stop("Please set the number of factors before attempting to pass additional arguments along to fsvsim.")
  }

  if ("a" %in% prior_obj$supplied_args) {
    if (!is.atomic(prior_obj$a) || length(prior_obj$a) > 1) {
      stop("a must be a vector with a single element.")
    }
  } else {
    prior_obj$supplied_args <- c(prior_obj$supplied_args, "a")
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
#' @return  No return value, called for side effects.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' print(prior_obj)
print.mfbvar_prior <- function(x, ...) {
  cat("The following elements of the prior object have not been set: \n", names(sapply(x, is.null))[sapply(x, is.null)])
  cat("\n\n")
  cat("Checking if the steady-state prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n")
  }

  cat("Checking if a Minnesota-style prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
  }

  #cat("Checking if the Dirichlet-Laplace prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps) && !is.null(x$a)) {
    #cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps", "a")]
    #cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
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

  cat("\n")
}

#' Summary method for mfbvar_prior
#'
#' summary method for object of class mfbvar_prior, showing some basic
#' information regarding the contents of the prior.
#' @param object prior object (class \code{mfbvar_prior})
#' @param ... additional arguments (currently unused)
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{estimate_mfbvar}}, \code{\link{print.mfbvar_prior}}
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' summary(prior_obj)
summary.mfbvar_prior <- function(object, ...) {
  cat("PRIOR SUMMARY\n")
  cat("----------------------------\n")
  cat("General specification:\n")
  cat("  Y:", ncol(object$Y), "variables,", nrow(object$Y), "time points\n")
  cat("  aggregation:", object$aggregation, "\n")
  freq_count <- vapply(object$freqs, function(x, freq) sum(x == freq), numeric(1), freq = object$freq)
  freqs <- object$freqs
  freqs <- replace(freqs, freqs == "w", "weekly")
  freqs <- replace(freqs, freqs == "m", "monthly")
  freqs <- replace(freqs, freqs == "q", "quarterly")
  if (length(freq_count) == 1) {
    freq_cat <- sprintf("  freq: %d %s variables\n", freq_count, freqs)
  } else {
    freq_cat <- sprintf("  freq: %s variables\n", paste(sprintf("%d %s", rev(freq_count), rev(freqs)), collapse = ", "))
  }
  cat(freq_cat)
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
  cat("  block_exo:", ifelse(is.null(object$block_exo), "0", length(object$block_exo)), "block exogenous variables\n")
  cat("  n_lags:", object$n_lags, "\n")
  cat("  n_fcst:", object$n_fcst, "\n")
  cat("  n_burnin:", object$n_burnin, "\n")
  cat("  n_reps:", object$n_reps, "\n")
  cat("----------------------------\n")
  cat("Steady-state prior:\n")
  cat("  d:", ifelse(is.null(object$d), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(ncol(object$d), "deterministic variables"))),"\n")
  cat("  d_fcst:", ifelse(is.null(object$d_fcst), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(nrow(object$d_fcst), "forecasts, ", ncol(object$d), "deterministic variables"))),"\n")
  cat("  prior_psi_mean:", ifelse(is.null(object$prior_psi_mean), "<missing>", "prior mean vector of steady states"), "\n")
  cat("  prior_psi_Omega:", ifelse(is.null(object$prior_psi_Omega), "<missing>", "prior covariance matrix of steady states"), "\n")
  cat("  check_roots:", object$check_roots, "\n")
  cat("----------------------------\n")
  cat("Hierarchical steady-state prior:\n")
  cat("  s:", ifelse(is.null(object$s), "<missing>", object$s), "\n")
  cat("  c0:", ifelse(is.null(object$prior_ng), "<missing>", object$prior_ng[1]), "\n")
  cat("  c1:", ifelse(is.null(object$prior_ng), "<missing>", object$prior_ng[2]), "\n")
  cat("----------------------------\n")
  #cat("Dirichlet-Laplace prior:\n")
  #cat("  a:", ifelse(is.null(object[["a"]]), "<missing>", object[["a"]]), "\n")
  #cat("----------------------------\n")
  cat("Common stochastic volatility:\n")
  cat(sprintf("  prior_phi: mean = %g, var = %g", object$prior_phi[1], object$prior_phi[2]), "\n")
  cat(sprintf("  prior_sigma2: mean = %g, df = %d", object$prior_sigma2[1], object$prior_sigma2[2]), "\n")
  cat("----------------------------\n")
  cat("Factor stochastic volatility:\n")
  cat("  n_fac:", ifelse(is.null(object$n_fac), "<missing>", object$n_fac), "\n")
  cat("  n_cores:", ifelse(is.null(object$n_cores), "<missing>", object$n_cores), "\n")
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
  if ("restrict" %in% object$supplied_args) {
    cat("  restrict:", object$restrict, "\n")
  }
  cat("----------------------------\n")
  cat("Other:\n")
  cat("  verbose:", object$verbose, "\n")

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

  if (!(prior %in% c("ss", "ssng", "minn", "dl"))) {
    stop("prior must be 'ss', 'ssng', 'minn' or 'dl'.")
  }
  if (!(variance %in% c("iw", "fsv", "csv", "diffuse"))) {
    stop("volatility must be 'iw', 'diffuse', 'csv' or 'fsv'.")
  }

  if (prior == "dl" && !(variance %in% c("fsv", "diffuse"))) {
    stop("The Dirichlet-Laplace prior (dl) can only be used with variance specifications fsv and diffuse.")
  }

  class(mfbvar_prior) <- c(sprintf("mfbvar_%s_%s", prior, variance), sprintf("mfbvar_%s", prior), sprintf("mfbvar_%s", variance), class(mfbvar_prior))

  time_out <- c(time_out, Sys.time())
  main_run <-  mcmc_sampler(mfbvar_prior)
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
    dimnames(main_run$psi) <- list(iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin),
                                   param = paste0(rep(names_col, n_determ), ".", rep(names_determ, each = n_vars)))
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars)),
                                  iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
  } else {
    dimnames(main_run$Pi) <- list(dep = names_col,
                                  indep = c("const", paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars))),
                                  iteration = 1:(mfbvar_prior$n_reps/mfbvar_prior$n_thin))
  }

  if (sum(mfbvar_prior$freq == "m") == 0 || sum(mfbvar_prior$freq == "m") == ncol(mfbvar_prior$Y)) {
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


#' Printing method for class mfbvar
#'
#' Method for printing \code{mfbvar} objects.
#'
#' @param x object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @template man_template
#' @return No return value, called for side effects.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' mod_minn

print.mfbvar <- function(x, ...){
  ss <- ifelse(x$prior == "ss", "steady-state ", "")
  freq_type <- ifelse(sum(x$freq == "m") == 0, "Quarterly", ifelse(sum(x$freq == "q") == 0, "Monthly", "Mixed-frequency"))
  var_type <- switch(x$variance,
                iw = "Inverse Wishart",
                diffuse = "Diffuse",
                fsv = sprintf("Factor stochastic volatility (%d factors)", x$mfbvar_prior$n_fac),
                csv = "Common stochastic volatility")
  cat(paste0(sprintf("%s BVAR with:\n", freq_type), ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "),
             "\nPrior: ", x$prior, "\n",
             "\nError covariance matrix: ", var_type, "\n",
             x$n_lags, " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecast\n",
             x$n_reps, " draws used in main chain"))
  cat("\n")
}

#' Summary method for class mfbvar
#'
#' Method for summarizing \code{mfbvar} objects.
#'
#' @param object object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @template man_template
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' summary(mod_minn)

summary.mfbvar <- function(object, ...){
  print(object)
}

#' Plotting methods for posterior mfbvar objects
#'
#' Methods for plotting posterior mfbvar objects.
#' @param x object of class \code{mfbvar_minn} or \code{mfbvar_ss}
#' @param aggregate_fcst Boolean indicating whether forecasts of the latent monthly series should be aggregated to the quarterly frequency.
#' @param plot_start Time period (date or number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param variables Vector of names or positions of variables to include in the plot of variances
#' @param pred_bands Single number (between \code{0.0} and \code{1.0}) giving the coverage level of forecast intervals.
#' @param ss_bands (Steady-state prior only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior steady-state intervals.
#' @param var_bands (\code{varplot} only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior intervals for the error standard deviations.
#' @param nrow_facet an integer giving the number of rows to use in the facet
#' @param ... Currently not in use.
#' @return A \code{\link[ggplot2]{ggplot}}.
#' @name plot-mfbvar
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, d = "intercept",
#'                        n_lags = 4, n_reps = 20,
#'                        n_fcst = 4, n_fac = 1)
#'
#' prior_intervals <- matrix(c(1, 3,
#'                             4, 8,
#'                             1, 3), ncol = 2, byrow = TRUE)
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
plot.mfbvar_ss <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ss_bands = 0.95, ...){


  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    stop("To plot the forecasts, proper dates must be provided in the input data.")
  }
  fcst_start <- lubridate::as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)


  plot_range_names <- fcst_start %m+% months(-nrow(x$Y):(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(nrow(x$Y)-x$n_fcst*5, 0):nrow(x$Y)
    }  else {
      plot_range <- 1:nrow(x$Y)
    }
  } else {
    plot_start <- tryCatch(lubridate::as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):nrow(x$Y)
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
  names_row <- if (is.null(x$names_row)) 1:nrow(x$Y) else x$names_row
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
    preds <- predict(x, aggregate_fcst = aggregate_fcst, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(data.frame(variable = names_col[i],
                                     time = (1:nrow(x$Y))[last_pos[i]],
                                     fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]),
                    fcst)
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
  breaks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  breaks <- na.omit(breaks)
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
plot.mfbvar_ssng <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ss_bands = 0.95, ...) {
  plot.mfbvar_ss(x, aggregate_fcst = aggregate_fcst, plot_start = plot_start,
                 pred_bands = pred_bands, nrow_facet = nrow_facet, ss_bands = ss_bands, ...)
}

#' @rdname plot-mfbvar
plot.mfbvar_minn <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                             pred_bands = 0.8, nrow_facet = NULL, ...){

  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    stop("To plot the forecasts, proper dates must be provided in the input data.")
  }
  fcst_start <-lubridate::as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)


  plot_range_names <- fcst_start %m+% months(-nrow(x$Y):(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(nrow(x$Y)-x$n_fcst*5, 0):nrow(x$Y)
    }  else {
      plot_range <- 1:nrow(x$Y)
    }
  } else {
    plot_start <- tryCatch(lubridate::as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):nrow(x$Y)
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
    preds <- predict(x, aggregate_fcst = aggregate_fcst, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(variable = names_col[i],
                                     time = (1:nrow(x$Y))[last_pos[i]],
                                     fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
                                     lower = x$Y[last_pos[i], i],
                                     median = x$Y[last_pos[i], i],
                                     upper  = x$Y[last_pos[i], i]),
                    fcst)
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
  breaks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  breaks <- na.omit(breaks)
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

plot.mfbvar_dl <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ...) {
  plot.mfbvar_minn(x, aggregate_fcst = aggregate_fcst, plot_start = plot_start,
                   pred_bands = pred_bands, nrow_facet = nrow_facet, ...)
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
    variances_fsv(variances, x$h, x$facload, variables_num, n_fac, n_reps, n_T, n_vars, n_plotvars)
  }
  if (sv_type == "csv") {
    variances_csv(variances, x$Sigma, x$f, n_T, n_reps, variables_num)
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





#' Plot method for class \code{mfbvar_prior}
#'
#' Method for plotting \code{mfbvar_prior} objects.
#'
#' @param x object of class \code{mfbvar_prior}
#' @param nrow_facet number of rows in facet
#' @param ... Currently not in use.
#' @details The function plots the data. If the prior moments for the steady-state parameters are available in \code{x}, these are included.
#' @return A \code{\link[ggplot2]{ggplot}}.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20, n_fcst = 4)
#' plot(prior_obj)
plot.mfbvar_prior <- function(x, nrow_facet = NULL, ...){


  ss_level <- c(0.025, 0.975)

  names_col <- if (is.null(colnames(x$Y))) paste0("x", 1:x$n_vars) else colnames(x$Y)

  if (!is.null(x$d)& !is.null(x$prior_psi_mean) & !is.null(x$prior_psi_Omega)) {
    ss_flag <- TRUE
  } else {
    ss_flag <- FALSE
  }

  if (!is.null(x$d)& !is.null(x$prior_psi_mean)) {
    ssng_flag <- TRUE
  } else {
    ssng_flag <- FALSE
  }

  if (ss_flag) {
    n_determ <- ncol(x$d)
    ss_lower  <- x$d %*% t(matrix(qnorm(ss_level[1], x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))
    ss_median <- x$d %*% t(matrix(qnorm(0.5, x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))
    ss_upper  <- x$d %*% t(matrix(qnorm(ss_level[2], x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))

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

