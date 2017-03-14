#' Interval to moments
#'
#' Convert a matrix of 95 \% prior probability intervals for the steady states to prior moments.
#' @templateVar prior_psi_int TRUE
#' @template man_template
#' @return A list with two components:
#' \item{prior_psi_mean}{The prior mean of psi}
#' \item{prior_psi_Omega}{The prior covariance matrix of psi}

interval_to_moments <- function(prior_psi_int) {
  stopifnot(is.matrix(prior_psi_int), ncol(prior_psi_int) == 2)
  prior_psi_mean <- rowMeans(prior_psi_int)
  prior_psi_Omega <- diag(((prior_psi_int[, 2] - prior_psi_int[, 1]) / (qnorm(0.975, mean = 0, sd = 1)*2))^2)
  return(list(prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega))
}
