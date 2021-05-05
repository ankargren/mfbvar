#' Interval to moments
#'
#' Convert a matrix of \code{100*(1-alpha)} \% prior probability intervals for the steady states to prior moments.
#' @param prior_psi_int Matrix of size \code{(n_determ * n_vars) x 2} with the prior intervals
#' @param alpha \code{100*(1-alpha)} is the prior probability of the interval
#' @return A list with two components:
#' \item{prior_psi_mean}{The prior mean of psi}
#' \item{prior_psi_Omega}{The prior covariance matrix of psi}
#' @examples
#' prior_intervals <- matrix(c(
#'   0.1, 0.2,
#'   0.4, 0.6
#' ), ncol = 2, byrow = TRUE)
#' psi_moments <- interval_to_moments(prior_intervals)
interval_to_moments <- function(prior_psi_int, alpha = 0.05) {
  stopifnot(is.matrix(prior_psi_int), ncol(prior_psi_int) == 2)
  prior_psi_mean <- rowMeans(prior_psi_int)
  prior_psi_Omega <- diag(((prior_psi_int[, 2] - prior_psi_int[, 1]) / (qnorm(1 - alpha / 2, mean = 0, sd = 1) * 2))^2)
  return(list(prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega))
}
