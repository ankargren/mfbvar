#' @title OLS functions
#' @description Helper functions for multivariate regression and sum of squared error computations
#' @param X The regressor matrix.
#' @param Y The dependnet variable matrix.
#' @keywords internal
#' @noRd
#' @return
#' \item{pi_sample}{Estimated coefficients.}
ols_pi <- function(X, Y) {
  ridge <- 1e-6
  error_count <- 0
  fail <- TRUE
  while (fail) {
    pi_sample <- tryCatch(
      {
        solve(crossprod(X) + diag(ridge, ncol(X))) %*% crossprod(X, Y)
      },
      error = function(cond) {
        cond
      }
    )
    if (!inherits(pi_sample, "error")) {
      fail <- FALSE
    } else {
      ridge <- ridge * 10
    }
  }
  return(pi_sample)
}

#' @rdname ols_pi
#' @param Pi The estimated coefficients.
#' @keywords internal
#' @noRd
#' @return
#' \item{s_sample}{The sum of squared residuals matrix.}
ols_s <- function(X, Y, Pi) {
  s_sample <- crossprod(Y - X %*% Pi)
  return(s_sample)
}

#' Initialize Gibbs sampler using OLS
#'
#' Initializes the Gibbs sampler using OLS.
#' @param z A matrix of size \code{(n_T + n_lags) x n_vars} of data
#' @param d The matrix of size \code{(n_T + n_lags) x n_determ} of deterministic terms
#' @param n_lags Number of lags
#' @param n_T Number of time periods
#' @param Number of variables
#' @param Number of deterministic components
#' @return A list with components:
#' \item{Gam}{A matrix of size \code{n_vars x (n_vars * n_lags + n_determ)} of estimated parameters.}
#' \item{S}{Estimated error covariance matrix.}
#' \item{psi}{The estimated steady-state parameters.}
#' @keywords internal
#' @noRd

ols_initialization <- function(z, d, n_lags, n_T, n_vars, n_determ) {
  n_T <- nrow(z)
  # Create regressor matrix (this is z in Karlsson, 2013)
  XX <- c()
  for (i in 1:n_lags) {
    XX <- cbind(XX, z[(n_lags + 1 - i):(n_T - i), ])
  }
  XX <- cbind(XX, d[(n_lags + 1):n_T, ])
  YY <- z[(n_lags + 1):n_T, ]

  # Gamma in Karlsson (2013, p. 797)
  Gam <- t(ols_pi(XX, YY))
  Pi <- Gam[, 1:(n_vars * n_lags)]
  const <- Gam[, (n_vars * n_lags + 1):(n_vars * n_lags + n_determ)]
  psi <- c(solve(diag(n_vars) - Pi %*%
    kronecker(matrix(1, n_lags, 1), diag(n_vars))) %*% const)

  return(list(
    Pi = Pi, S = crossprod(YY - XX %*% t(Gam)) / n_T,
    psi = psi, const = const
  ))
}
