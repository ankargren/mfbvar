#' OLS functions
#'
#' Helper functions for multivariate regression and sum of squared error computations
#'
#' @param X The regressor matrix.
#' @param Y The dependnet variable matrix.
#' @return
#' \item{pi_sample}{Estimated coefficients.}

ols_pi <- function(X, Y) {
  pi_sample <- chol2inv(chol(crossprod(X))) %*% crossprod(X, Y)
  return(pi_sample)
}

#' @rdname ols_pi
#' @param Pi The estimated coefficients.
#' @return
#' \item{s_sample}{The sum of squared residuals matrix.}
ols_s <- function(X, Y, Pi) {
  s_sample <- crossprod(Y - X %*% Pi)
  return(s_sample)
}

#' Initialize Gibbs sampler using OLS
#'
#' Initializes the Gibbs sampler using OLS.
#' @templateVar z TRUE
#' @templateVar d TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_T TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_determ TRUE
#' @template man_template
#' @return A list with components:
#' \item{Gam}{A matrix of size \code{n_vars * (n_vars*n_lags +n_determ)} of estimated parameters.}
#' \item{S}{Estimated error covariance matrix.}
#' \item{psi}{The estimated steady-state parameters.}

ols_initialization <- function(z, d, n_lags, n_T, n_vars, n_determ) {
  n_T <- nrow(z)
  # Create regressor matrix (this is z in Karlsson, 2013)
  XX <- c()
  for (i in 1:n_lags) {
    XX <- cbind(XX, z[(n_lags+1-i):(n_T - i), ])
  }
  XX <- cbind(XX, d[(n_lags+1):n_T, ])
  YY <- z[(n_lags+1):n_T, ]

  # Gamma in Karlsson (2013, p. 797)
  Gam <- t(ols_pi(XX, YY))
  Pi  <- Gam[, 1:(n_vars * n_lags)]
  psi <- c(solve(diag(n_vars) - Pi %*%
                   kronecker(matrix(1, n_lags, 1), diag(n_vars))) %*%
             Gam[, (n_vars * n_lags + 1):(n_vars * n_lags + n_determ)])
  return(list(Pi = Pi, S = crossprod(YY - XX %*% t(Gam)) / n_T,
              psi = psi))
}
