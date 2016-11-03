# Contains:
# posterior_psi
# posterior_psi_omega
# posterior_s
# posterior_pi
# posterior_pi_omega

posterior_psi <- function(psi_omega, U, sigma, Y_tilde, D, prior_psi_omega, prior_psi) {
  sigmaYD <- solve(sigma) %*% t(Y_tilde) %*% D
  dim(sigmaYD) <- c(prod(dim(sigmaYD)), 1)
  psi <- psi_omega %*% (t(U) %*% sigmaYD + solve(prior_psi_omega) %*% prior_psi)
  return(psi)
}

posterior_psi_omega <- function(U, D, sigma, prior_psi_omega) {
  psi_omega <- solve(t(U) %*% (kronecker(crossprod(D), solve(sigma))) %*% U + solve(prior_psi_omega))
  return(psi_omega)
}

posterior_s <- function(prior_s, s_sample, prior_pi, pi_sample, prior_pi_omega, Z) {
  s <- prior_s + s_sample + t(prior_pi - pi_sample) %*% solve(prior_pi_omega + solve(crossprod(Z))) %*% (prior_pi - pi_sample)
}

posterior_pi <- function(pi_omega, prior_pi_omega, prior_pi, demeaned_Z) {
  n_vars <- dim(prior_pi)[2]
  Z_tilde_1T1 <- demeaned_Z[-nrow(demeaned_Z), ]
  z_tilde_2T  <- demeaned_Z[-1, 1:n_vars]
  pi <- pi_omega %*% (solve(prior_pi_omega) %*% prior_pi + crossprod(Z_tilde_1T1, z_tilde_2T))
  return(pi)
}

posterior_pi_omega <- function(prior_pi_omega, Z_tilde_1T1) {
  pi_omega <- solve(solve(prior_pi_omega) + crossprod(Z_tilde_1T1))
  return(pi_omega)
}
