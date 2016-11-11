# Contains:
# posterior_psi
# posterior_psi_omega
# posterior_s
# posterior_pi
# posterior_pi_omega

posterior_psi <- function(U, D, sigma, prior_psi_omega, psi_omega, Y_tilde, prior_psi) {
  sigmaYD <- matrix(c(solve(sigma) %*% t(Y_tilde) %*% D), ncol = 1)
  psi <- psi_omega %*% (t(U) %*% sigmaYD + solve(prior_psi_omega) %*% prior_psi)
  return(psi)
}

posterior_psi_omega <- function(U, D, sigma, prior_psi_omega) {
  psi_omega <- solve(t(U) %*% (kronecker(crossprod(D), solve(sigma))) %*% U + solve(prior_psi_omega))
  return(psi_omega)
}
