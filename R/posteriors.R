#' Draw from posterior of Pi and Sigma
#'
#' Function for drawing from the posterior of Pi and Sigma, which can be used as a block in a Gibbs sampler.
#' @templateVar Z_r1 TRUE
#' @templateVar d TRUE
#' @templateVar psi_r1 TRUE
#' @templateVar prior_Pi_mean TRUE
#' @templateVar prior_Pi_Omega TRUE
#' @templateVar inv_prior_Pi_Omega TRUE
#' @templateVar Omega_Pi TRUE
#' @templateVar prior_S TRUE
#' @templateVar prior_nu TRUE
#' @templateVar check_roots TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_T TRUE
#' @template man_template
#' @keywords internal
#' @return \code{posterior_Pi_Sigma} returns a list with:
#' \item{Pi_r}{The draw of \code{Pi}.}
#' \item{Sigma_r}{The draw of \code{Sigma}.}
#' \item{num_try}{The try at which a stable draw was obtained.}
#' \item{root}{The maximum eigenvalue (in modulus) of the system.}
posterior_Pi_Sigma <- function(Z_r1, d, psi_r1, prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T) {
  ################################################################
  ### Preliminary calculations

  # Demean z, create Z (companion form version)
  demeaned_z <- Z_r1 - d %*% t(matrix(psi_r1, nrow = n_vars))
  demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
  XX <- demeaned_Z[-nrow(demeaned_Z), ]
  YY <- demeaned_Z[-1, 1:n_vars]
  XXt.XX <- crossprod(XX)
  XXt.XX.inv <- chol2inv(chol(XXt.XX))
  Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)

  ################################################################
  ### Pi and Sigma step

  # Posterior moments of Pi
  post_Pi_Omega <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
  post_Pi       <- post_Pi_Omega %*% (Omega_Pi + crossprod(XX, YY))

  # Then Sigma
  s_sample  <- crossprod(YY - XX %*% Pi_sample)
  Pi_diff <- prior_Pi_mean - Pi_sample
  post_s <- prior_S + s_sample + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff
  post_nu <- n_T + prior_nu
  Sigma_r <- rinvwish(v = post_nu, S = post_s)


  # Draw Pi conditional on Sigma
  # This ensures that the draw is stationary
  stationarity_check <- FALSE
  iter <- 0
  Pi_temp <- array(NA, dim = c(n_vars, n_vars * n_lags, ifelse(check_roots, 1000, 1)))
  while(stationarity_check == FALSE) {
    iter <- iter + 1
    Pi_temp[,,iter] <- rmatn(M = t(post_Pi), Q = post_Pi_Omega, P = Sigma_r)
    Pi_comp    <- build_companion(Pi_temp[,, iter], n_vars = n_vars, n_lags = n_lags)
    if (check_roots == TRUE) {
      root <- max_eig_cpp(Pi_comp)
    } else {
      root <- 0
    }
    if (root < 1) {
      stationarity_check <- TRUE
      num_try <- iter
      Pi_r <- Pi_temp[,,iter]
    }
    if (iter == 1000) {
      stop("Attempted to draw stationary Pi 1,000 times.")
    }
  }

  return(list(Pi_r = Pi_r, Sigma_r = Sigma_r, num_try = num_try, root = root))
}

#' Draw from posterior of psi
#'
#' Function for drawing from the posterior of psi, which can be used as a block in a Gibbs sampler.
#' @inherit posterior_Pi_Sigma
#' @templateVar Pi_r TRUE
#' @templateVar Sigma_r TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar D_mat TRUE
#' @templateVar n_determ TRUE
#' @template man_template
#' @keywords internal
#' @return \code{posterior_psi} returns:
#' \item{psi_r}{The draw of \code{psi}.}
posterior_psi <- function(Pi_r, Sigma_r, Z_r1, prior_psi_mean, prior_psi_Omega, D_mat, n_vars, n_lags, n_determ) {
  U <- build_U_cpp(Pi = Pi_r, n_determ = n_determ,
                   n_vars = n_vars, n_lags = n_lags)
  post_psi_Omega <- posterior_psi_Omega(U = U, D_mat = D_mat, Sigma = Sigma_r,
                                        prior_psi_Omega = prior_psi_Omega)
  Y_tilde <- build_Y_tilde(Pi = Pi_r, z = Z_r1)

  post_psi <- posterior_psi_mean(U = U, D_mat = D_mat, Sigma = Sigma_r, prior_psi_Omega = prior_psi_Omega,
                                 post_psi_Omega = post_psi_Omega, Y_tilde = Y_tilde, prior_psi_mean = prior_psi_mean)
  psi_r <- t(rmultn(m = post_psi, Sigma = post_psi_Omega))
  return(psi_r)
}

#' Compute posterior moments of the steady-state parameters
#'
#' Computes the mean and variance of the conditional posterior distribution of the steady-state parameters.
#' @templateVar U TRUE
#' @templateVar D_mat TRUE
#' @templateVar Sigma TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar post_psi_Omega TRUE
#' @templateVar Y_tilde TRUE
#' @templateVar prior_psi_mean TRUE
#' @template man_template
#' @keywords internal
#' @return The return is:
#' \item{psi}{The posterior mean (from \code{\link{posterior_psi_mean}})}
posterior_psi_mean <- function(U, D_mat, Sigma, prior_psi_Omega, post_psi_Omega, Y_tilde, prior_psi_mean) {
  SigmaYD <- matrix(c(chol2inv(chol(Sigma)) %*% t(Y_tilde) %*% D_mat), ncol = 1)
  psi <- post_psi_Omega %*% (t(U) %*% SigmaYD + chol2inv(chol(prior_psi_Omega)) %*% prior_psi_mean)
  return(psi)
}

#' @rdname posterior_psi_mean
#' @keywords internal
#' @return \item{psi_Omega}{The posterior variance (from \code{\link{posterior_psi_Omega}})}
posterior_psi_Omega <- function(U, D_mat, Sigma, prior_psi_Omega) {
  psi_Omega <- chol2inv(chol(t(U) %*% (kronecker(crossprod(D_mat), chol2inv(chol(Sigma)))) %*% U + chol2inv(chol(prior_psi_Omega))))
  return(psi_Omega)
}

