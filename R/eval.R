#' Evaluate the conditional posterior of Pi and Sigma using Rao-Blackwellization
#'
#' Evaluates the conditional posterior of Pi and Sigma using Rao-Blackwellization of the draws from the Gibbs sampler.
#' @templateVar Z_array TRUE
#' @templateVar d TRUE
#' @templateVar post_psi_center TRUE
#' @templateVar post_Pi_center TRUE
#' @templateVar post_Sigma_center TRUE
#' @templateVar post_nu TRUE
#' @templateVar prior_Pi TRUE
#' @templateVar prior_Pi_Omega TRUE
#' @templateVar prior_S TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUEF
#' @templateVar n_reps TRUE
#' @template man_template
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
    demeaned_z <- Z_array[,,i+1] - d %*% post_psi_center
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    XXt.XX <- crossprod(XX)
    XXt.XX.inv <- chol2inv(chol(XXt.XX))
    Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_Pi_Omega_i <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
    post_Pi_i       <- post_Pi_Omega_i %*% (Omega_Pi + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% Pi_sample)
    Pi_diff <- prior_Pi_mean - Pi_sample
    post_s_i <- prior_S + s_sample + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff

    # Evaluate
    evals[i] <- exp(dnorminvwish(X = t(post_Pi_center), Sigma = post_Sigma_center, M = post_Pi_i, P = post_Pi_Omega_i, S = post_s_i, v = post_nu))
  }

  return(evals)

}

#' Evaluate the marginal posterior of psi
#'
#' Evaluates the marginal posterior of psi using the draws from the Gibbs sampler.
#' @templateVar Pi_array TRUE
#' @templateVar Sigma_array TRUE
#' @templateVar Z_array TRUE
#' @templateVar post_psi_center TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar D_mat TRUE
#' @templateVar n_determ TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_reps TRUE
#' @template man_template
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
    U <- build_U_cpp(Pi = Pi_array[,,r], n_determ = n_determ,
                     n_vars = n_vars, n_lags = n_lags)
    post_psi_Omega <- posterior_psi_Omega(U = U, D_mat = D_mat, Sigma = Sigma_array[,, r],
                                          prior_psi_Omega = prior_psi_Omega)
    Y_tilde <- build_Y_tilde(Pi = Pi_array[,, r], z = Z_array[,, r])

    post_psi_center <- posterior_psi_mean(U = U, D_mat = D_mat, Sigma = Sigma_array[,, r], prior_psi_Omega = prior_psi_Omega,
                                          post_psi_Omega = post_psi_Omega, Y_tilde = Y_tilde, prior_psi_mean = prior_psi_mean)

    evals[r] <- exp(dmultn(x = post_psi_center, m = post_psi_center, Sigma = post_psi_Omega))
  }

  return(evals)
}
