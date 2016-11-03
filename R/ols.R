ols_pi <- function(X, Y) {
  pi_sample <- solve(crossprod(X)) %*% crossprod(X, Y)
  return(pi_sample)
}

ols_s <- function(X, Y, Pi) {
  s_sample <- crossprod(Y - X %*% Pi)
  return(s_sample)
}


ols_initialization <- function(z, d, n_lags) {
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
  return(list(Gam = Gam, S = crossprod(YY - XX %*% t(Gam)) / n_T))
}
