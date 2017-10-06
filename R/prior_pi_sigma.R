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

  prior_S <- 1/(prior_nu - n_vars - 1) * diag(error_variance)


  return(list(prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = diag(prior_Pi_Omega), prior_S = prior_S))
}



