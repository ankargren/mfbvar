#' Create the priors for Pi and Sigma
#'
#' Creates the prior mean and covariance for Pi given the hyperparameters, and the prior parameters for Sigma.
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar Y TRUE
#' @templateVar n_lags TRUE
#' @templateVar prior_nu TRUE
#' @template man_template
#' @return \item{prior_Pi}{The prior mean matrix for Pi.}
#' \item{prior_Pi_Omega}{The prior covariance matrix for Pi.}
#' \item{prior_s}{The prior for Sigma.}
#' @keywords internal
#' @noRd
prior_Pi_Sigma <- function(lambda1, lambda2, prior_Pi_AR1, Y, n_lags, prior_nu) {
  # lambda1: 1-long vector (overall tightness)
  # lambda2: 1-long vector (lag decay)
  # prior_Pi_AR1: p-long vector with prior means for the AR(1) coefficients
  # Y: Txp matrix with data

  n_vars <- length(prior_Pi_AR1)
  prior_Pi_mean <- rbind(diag(prior_Pi_AR1), matrix(0, nrow = n_vars*(n_lags-1), ncol = n_vars))

  prior_Pi_Omega <- rep(0, n_lags * n_vars)
  error_variance <- rep(NA, n_vars)
  for (i in 1:n_vars) {
    success <- NULL
    init_order <- 4
    while(is.null(success)) {
      error_variance[i] <- tryCatch(arima(na.omit(Y[,i]), order = c(init_order, 0, 0), method = "ML")$sigma2,
                                 error = function(cond) NA)
      if (!is.na(error_variance[i])) {
        success <- 1
      } else {
        init_order <- init_order - 1
        if (init_order < 1) {
          stop("Too low order.")
        }
      }
    }
  }

  for (l in 1:n_lags) {
    for (r in 1:n_vars) {
      i <- (l - 1) * n_vars + r
      prior_Pi_Omega[i] <- lambda1^2 / (l^(lambda2) * sqrt(error_variance[r]))^2
    }
  }

  prior_S <- (prior_nu - n_vars - 1) * diag(error_variance)


  return(list(prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = diag(prior_Pi_Omega), prior_S = prior_S))
}

#' Create the priors for Pi and Sigma
#'
#' Creates the prior mean and covariance for Pi given the hyperparameters, and the prior parameters for Sigma.
#' @param lambda1 overall tightness
#' @param lambda2 cross-equation tightness
#' @param lambda3 lag decay
#' @param prior_Pi_AR1 prior means for AR(1) coefficients
#' @param Y data
#' @param n_lags number of lags
#' @return
#' \item{prior_Pi_Omega}{The prior covariance matrix for Pi.}
#' @keywords internal
#' @noRd
create_prior_Pi_Omega <- function(lambda1, lambda2, lambda3, prior_Pi_AR1, Y, n_lags,
                                  block_exo = NULL) {
  # lambda1: 1-long vector (overall tightness)
  # lambda2: 1-long vector (lag decay)
  # prior_Pi_AR1: p-long vector with prior means for the AR(1) coefficients
  # Y: Txp matrix with data

  n_vars <- length(prior_Pi_AR1)

  prior_Pi_Omega <- matrix(0, n_vars * n_lags + 1, n_vars)
  error_variance <- compute_error_variances(Y)

  prior_Pi_Omega[1, ] <- lambda1 * 100 * sqrt(error_variance)
  for (i in 1:n_vars) {
    prior_Pi_Omega[-1, i] <- lambda1 * lambda2^((1:n_vars) != i) * sqrt(error_variance[i])/
      ((rep(1:n_lags, each = n_vars))^(lambda3) * rep(sqrt(error_variance), times = n_lags))
    if (i %in% block_exo) {
      prior_Pi_Omega[-1, i] <- prior_Pi_Omega[-1, i] * (1e-06)^(!(1:ncol(Y) %in% block_exo))
    }
  }

  return(prior_Pi_Omega = prior_Pi_Omega^2)
}


