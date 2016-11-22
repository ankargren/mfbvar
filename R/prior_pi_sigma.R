prior_pi_sigma <- function(lambda1, lambda2, prior_mean, Y, n_lags, nu) {
  # lambda1: 1-long vector (overall tightness)
  # lambda2: 1-long vector (lag decay)
  # prior_mean: p-long vector with prior means for the AR(1) coefficients
  # Y: Txp matrix with data

  n_vars <- length(prior_mean)
  prior_pi <- rbind(diag(prior_mean), matrix(0, nrow = n_vars*(n_lags-1), ncol = n_vars))

  error_variance <- apply(Y, 2, function(x) arima(na.omit(x), order = c(4, 0, 0))$sigma^2)
  prior_pi_omega <- rep(0, n_lags * n_vars)
  for (l in 1:n_lags) {
    for (r in 1:n_vars) {
      i <- (l - 1) * n_vars + r
      prior_pi_omega[i] <- lambda1^2 / (l^(lambda2) * sqrt(error_variance[r]))^2
    }
  }

  prior_s <- 1/(nu - n_vars - 1) * diag(error_variance)


  return(list(prior_pi = prior_pi, prior_pi_omega = diag(prior_pi_omega), prior_s = prior_s))
}



