#' Create the priors for Pi and Sigma
#'
#' Creates the prior mean and covariance for Pi given the hyperparameters, and the prior parameters for Sigma.
#' @param lambda1 hyperparameter for overall tightness
#' @param lambda2 hyperparameter for cross-variable tightness (only if \code{independent = TRUE})
#' @param lambda3 hyperparameter for lag decay
#' @param lambda4 hyperparameter for intercept
#' @param prior_Pi_AR1 prior means for AR(1) coefficients
#' @param Y data matrix
#' @param n_lags number of lags
#' @param intercept boolean flag for indicating whether prior for intercept should be included
#' @param prior_nu prior degrees of freedom (only if \code{independent = FALSE})
#' @param block_exo vector of indexes of variables to be treated as block exogenous (only if \code{independent = TRUE})
#' @return The return object depends on whether the argument \code{independent = TRUE} or not. If \code{TRUE}, then the function returns a matrix of dimension (\code{n_vars * n_lags + intercept, n_vars}) in which each element is the prior variance of the corresponding element of the matrix of regression parameters Pi. If \code{independent = FALSE}, then the returned object is a list with components:
#' \describe{\item{prior_Pi_mean}{The prior mean matrix for Pi}
#' \item{prior_Pi_Omega}{The prior scale matrix for Pi.}
#' \item{prior_S}{The prior scale for Sigma}
#' \item{inv_prior_Pi_Omega}{Inverse of \code{prior_Pi_Omega}}
#' \item{Omega_Pi}{Precomputation of \code{inv_prior_Pi_Omega %*% prior_Pi_mean}}}
#' @keywords internal
#' @noRd
create_prior_Pi <- function(lambda1, lambda2,  lambda3, lambda4, prior_Pi_AR1, Y,
                           n_lags, intercept, prior_nu = NULL, block_exo = NULL,
                           independent) {
  n_vars <- length(prior_Pi_AR1)
  error_variance <- compute_error_variances(Y)

  prior_Pi_mean <- rbind(diag(prior_Pi_AR1), matrix(0, nrow = n_vars*(n_lags-1), ncol = n_vars))

  if (independent) {

    prior_Pi_Omega <- matrix(0, n_vars * n_lags + 1, n_vars)

    if (intercept) {
      prior_Pi_Omega[1, ] <- lambda1 * lambda4 * sqrt(error_variance)
      prior_Pi_mean <- rbind(0, prior_Pi_mean)
    }

    for (i in 1:n_vars) {
      prior_Pi_Omega[-1, i] <- lambda1 * lambda2^((1:n_vars) != i) * sqrt(error_variance[i])/
        ((rep(1:n_lags, each = n_vars))^(lambda3) * rep(sqrt(error_variance), times = n_lags))
      if (!is.null(block_exo) && (i %in% block_exo)) {
        prior_Pi_Omega[-1, i] <- prior_Pi_Omega[-1, i] * (1e-06)^(!(1:ncol(Y) %in% block_exo))
      }
    }

    if (!intercept) prior_Pi_Omega <- prior_Pi_Omega[-1, ]

    return(list(prior_Pi_Omega = prior_Pi_Omega^2,
                prior_Pi_mean = prior_Pi_mean))
  } else {

    prior_Pi_Omega <- rep(0, n_lags * n_vars)

    for (l in 1:n_lags) {
      for (r in 1:n_vars) {
        i <- (l - 1) * n_vars + r
        prior_Pi_Omega[i] <- lambda1^2 / (l^(lambda3) * sqrt(error_variance[r]))^2
      }
    }

    prior_S <- (prior_nu - n_vars - 1) * diag(error_variance)

    if (intercept) {
      # Add terms for constant
      prior_Pi_Omega <- c(lambda1^2*lambda4^2, prior_Pi_Omega)
      prior_Pi_mean <- rbind(0, prior_Pi_mean)
    }

    inv_prior_Pi_Omega <- chol2inv(chol(diag(prior_Pi_Omega)))
    Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

    return(list(prior_Pi_mean = prior_Pi_mean,
                prior_Pi_Omega = diag(prior_Pi_Omega),
                prior_S = prior_S,
                inv_prior_Pi_Omega = inv_prior_Pi_Omega,
                Omega_Pi = Omega_Pi))
  }
}

