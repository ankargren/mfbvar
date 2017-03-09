#' Create the priors for Pi and Sigma
#'
#' Creates the prior mean and covariance for Pi given the hyperparameters, and the prior parameters for Sigma.
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar Y TRUE
#' @templateVar n_lags TRUE
#' @templateVar nu TRUE
#' @template man_template
#' @return \item{prior_Pi}{The prior mean matrix for Pi.}
#' \item{prior_Pi_Omega}{The prior covariance matrix for Pi.}
#' \item{prior_s}{The prior for Sigma.}
prior_Pi_Sigma <- function(lambda1, lambda2, prior_Pi_AR1, Y, n_lags, nu) {
  # lambda1: 1-long vector (overall tightness)
  # lambda2: 1-long vector (lag decay)
  # prior_Pi_AR1: p-long vector with prior means for the AR(1) coefficients
  # Y: Txp matrix with data

  n_vars <- length(prior_Pi_AR1)
  prior_Pi_mean <- rbind(diag(prior_Pi_AR1), matrix(0, nrow = n_vars*(n_lags-1), ncol = n_vars))

  error_variance <- apply(Y, 2, function(x) arima(na.omit(x), order = c(4, 0, 0), method = "ML")$sigma2)
  prior_Pi_Omega <- rep(0, n_lags * n_vars)
  for (l in 1:n_lags) {
    for (r in 1:n_vars) {
      i <- (l - 1) * n_vars + r
      prior_Pi_Omega[i] <- lambda1^2 / (l^(lambda2) * sqrt(error_variance[r]))^2
    }
  }

  prior_s <- 1/(nu - n_vars - 1) * diag(error_variance)


  return(list(prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = diag(prior_Pi_Omega), prior_s = prior_s))
}



