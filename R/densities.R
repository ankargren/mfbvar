#' Normal inverse Wishart density function
#'
#' Density function for the (matrix) normal inverse Wishart distribution
#' @templateVar X TRUE
#' @templateVar Sigma TRUE
#' @templateVar M TRUE
#' @templateVar Q TRUE
#' @templateVar P TRUE
#' @templateVar S TRUE
#' @templateVar v TRUE
#' @keywords internal
#' @template man_template
#' @return
#' For \code{dnorminvwish}: the evaluated density.\\n
#' For \code{rmatn} or \code{rinvwish}: the random numbers.
dnorminvwish <- function(X, Sigma, M, P, S, v) {
  q <- dim(Sigma)[1]
  p <- dim(P)[1]
  det_Sigma <- det(Sigma)
  inv_Sigma <- chol2inv(chol(Sigma))
  dmultnorm <- (-p*q/2) * log(2 * pi) + (-p/2) * log(det_Sigma) + (-q/2)*log(det(P)) + (-1/2 * sum(diag(inv_Sigma %*% t(X - M) %*% chol2inv(chol(P)) %*% (X - M))))
  cc <- (v * q/2)*log(2) + (q*(q-1)/4)*log(pi) + sum(lgamma((v+1-1:q)/2))
  dinvwish <- -cc + (v/2) * log(det(S)) -(v+q+1)/2*log(det_Sigma) -1/2 * sum(diag(inv_Sigma %*% S))
  return(dmultnorm + dinvwish)
}

#' Multivariate normal density function
#'
#' Density function for the multivariate normal distribution
#' @templateVar x TRUE
#' @templateVar m TRUE
#' @template man_template
#' @inherit dnorminvwish
#' @keywords internal
#' @return
#' For \code{dmultn}: the evaluated density.\\n
#' For \code{rmultn}: \eqn{p} random numbers.
dmultn <- function(x, m, Sigma) {
  log_d <- (-1/2)* log(det(2*pi*Sigma)) -1/2 * t(x-m) %*% chol2inv(chol(Sigma)) %*% (x-m)
  return(log_d)
}

#' Truncated multivariate normal density function
#'
#' Density function for the truncated multivariate normal distribution
#' @templateVar V_inv TRUE
#' @param d The number of components.
#' @templateVar p_trunc TRUE
#' @templateVar chisq_val TRUE
#' @template man_template
#' @keywords internal
#' @inherit dmultn
dnorm_trunc <- function(x, m, V_inv, d, p_trunc, chisq_val) {
  qf <- t(x - m) %*% V_inv %*% (x - m)
  return((1/p_trunc) * (1/sqrt((2*pi)^d/det(V_inv))) * exp(-0.5 * qf) * (qf < chisq_val))
}

