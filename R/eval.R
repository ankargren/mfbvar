#' Evaluate the conditional posterior of Pi and Sigma using Rao-Blackwellization
#'
#' Evaluates the conditional posterior of Pi and Sigma using Rao-Blackwellization of the draws from the Gibbs sampler.
#' @param Z_array The array of draws of Z from the Gibbs sampler
#' @param d The matrix of size \code{(n_T + n_lags) * n_determ} of deterministic
#' terms
#' @param post_psi_center The value at which to do the evaluation (e.g. the
#' posterior mean/median)
#' @param post_Pi_center The value at which to do the evaluation (e.g. the
#' posterior mean/median)
#' @param post_Sigma_center The value at which to do the evaluation (e.g. the
#' posterior mean/median)
#' @param post_nu The posterior of the parameter \eqn{\nu}
#' @param prior_Pi Matrix of size \code{n_vars x (n_vars * n_lags)} containing
#' the prior for the mean of the dynamic coefficients.
#' @param prior_Pi_Omega Matrix of size \code{(n_vars * n_lags) x
#' (n_vars * n_lags)} containing the prior for (part of) the prior covariance of
#' the dynamic coefficients
#' @param prior_S The prior for \eqn{\Sigma}
#' @param n_vars Number of variables
#' @param n_lags Number of lags
#' @param n_reps Number of draws
#' @keywords internal
#' @noRd
#' @return The return is:
#' \item{evals}{A vector with the evaulations.}
#'
eval_Pi_Sigma_RaoBlack <- function(Z_array, d, post_psi_center, post_Pi_center, post_Sigma_center, post_nu, prior_Pi_mean, prior_Pi_Omega, prior_S, n_vars, n_lags, n_reps) {
  ################################################################
  ### Compute the Rao-Blackwellized estimate of posterior

  evals <- vector("numeric", n_reps - 1)
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean
  for (i in 1:length(evals)) {
    # Demean z, create Z_array (companion form version)
    demeaned_z <- Z_array[, , i + 1] - d %*% post_psi_center
    demeaned_Z <- mfbvar:::build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    XXt.XX <- crossprod(XX)
    XXt.XX.inv <- chol2inv(chol(XXt.XX))
    Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_Pi_Omega_i <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
    post_Pi_i <- post_Pi_Omega_i %*% (Omega_Pi + crossprod(XX, YY))

    # Then Sigma
    s_sample <- crossprod(YY - XX %*% Pi_sample)
    Pi_diff <- prior_Pi_mean - Pi_sample
    post_s_i <- prior_S + s_sample + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff

    # Evaluate
    evals[i] <- dnorminvwish(X = t(post_Pi_center), Sigma = post_Sigma_center, M = post_Pi_i, P = post_Pi_Omega_i, S = post_s_i, v = post_nu)
  }

  return(evals)
}

#' Evaluate the marginal posterior of psi
#'
#' Evaluates the marginal posterior of psi using the draws from the Gibbs sampler.
#' @param Pi_array Array of draws of Pi from the Gibbs sampler
#' @param Sigma_array Array of draws of Sigma from the Gibbs sampler
#' @param Z_array The array of draws of Z from the Gibbs sampler
#' @param post_psi_center The value at which to do the evaluation (e.g. the
#' posterior mean/median)
#' @param prior_psi_mean Vector of length \code{n_determ * n_vars} with the
#' prior means of the steady-state parameters
#' @param prior_psi_Omega Matrix of size \code{(n_determ * n_vars) x
#' (n_determ * n_vars)} with the prior covariance of the steady-state parameters
#' @param D_mat The \code{D} matrix (from \code{\link{build_DD}})
#' @param n_determ Number of deterministic components
#' @param n_vars Number of variables
#' @param n_lags Number of lags
#' @param n_reps Number of draws
#' @keywords internal
#' @noRd
#' @return The return is:
#' \item{evals}{A vector with the evaulations.}
#'
#'
eval_psi_MargPost <- function(Pi_array, Sigma_array, Z_array, post_psi_center, prior_psi_mean, prior_psi_Omega, D_mat, n_determ, n_vars, n_lags, n_reps) {
  post_psi_center <- matrix(post_psi_center, ncol = 1)

  evals <- vector("numeric", n_reps - 1)

  ################################################################
  ### Steady-state step
  for (r in 1:(n_reps - 1)) {
    U <- build_U_cpp(
      Pi = Pi_array[, , r], n_determ = n_determ,
      n_vars = n_vars, n_lags = n_lags
    )
    post_psi_Omega <- posterior_psi_Omega(
      U = U, D_mat = D_mat, Sigma = Sigma_array[, , r],
      prior_psi_Omega = prior_psi_Omega
    )
    Y_tilde <- build_Y_tilde(Pi = Pi_array[, , r], z = Z_array[, , r])

    post_psi_center <- posterior_psi_mean(
      U = U, D_mat = D_mat, Sigma = Sigma_array[, , r], prior_psi_Omega = prior_psi_Omega,
      post_psi_Omega = post_psi_Omega, Y_tilde = Y_tilde, prior_psi_mean = prior_psi_mean
    )

    evals[r] <- exp(dmultn(x = post_psi_center, m = post_psi_center, Sigma = post_psi_Omega))
  }

  return(evals)
}
