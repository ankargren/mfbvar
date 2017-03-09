#' Draw from posterior of Pi and Sigma
#'
#' Function for drawing from the posterior of Pi and Sigma, which can be used as a block in a Gibbs sampler.
#' @templateVar Z_r1 TRUE
#' @templateVar d TRUE
#' @templateVar psi_r1 TRUE
#' @templateVar prior_pi TRUE
#' @templateVar prior_pi_omega TRUE
#' @templateVar inv_prior_pi_omega TRUE
#' @templateVar omega_pi TRUE
#' @templateVar prior_s TRUE
#' @templateVar prior_nu TRUE
#' @templateVar check_roots TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_T TRUE
#' @template man_template
#' @return \code{pi_sigma_posterior} returns a list with:
#' \item{Pi_r}{The draw of \code{Pi}.}
#' \item{Sigma_r}{The draw of \code{Sigma}.}
#' \item{num_try}{The try at which a stable draw was obtained.}
#' \item{root}{The maximum eigenvalue (in modulus) of the system.}
pi_sigma_posterior <- function(Z_r1, d, psi_r1, prior_pi, prior_pi_omega, inv_prior_pi_omega, omega_pi, prior_s, prior_nu, check_roots, n_vars, n_lags, n_T) {
  ################################################################
  ### Preliminary calculations

  # Demean z, create Z (companion form version)
  demeaned_z <- Z_r1 - d %*% t(matrix(psi_r1, nrow = n_vars))
  demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
  XX <- demeaned_Z[-nrow(demeaned_Z), ]
  YY <- demeaned_Z[-1, 1:n_vars]
  XXt.XX <- crossprod(XX)
  XXt.XX.inv <- chol2inv(chol(XXt.XX))
  pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)

  ################################################################
  ### Pi and Sigma step

  # Posterior moments of Pi
  post_pi_omega <- solve(inv_prior_pi_omega + XXt.XX)
  post_pi       <- post_pi_omega %*% (omega_pi + crossprod(XX, YY))

  # Then Sigma
  s_sample  <- crossprod(YY - XX %*% pi_sample)
  pi_diff <- prior_pi - pi_sample
  post_s <- prior_s + s_sample + t(pi_diff) %*% solve(prior_pi_omega + XXt.XX.inv) %*% pi_diff
  nu <- n_T + prior_nu # Is this the right T? Or should it be T - lags?
  Sigma_r <- rinvwish(v = nu, S = post_s)


  # Draw Pi conditional on Sigma
  # This ensures that the draw is stationary
  stationarity_check <- FALSE
  iter <- 0
  Pi_temp <- array(NA, dim = c(n_vars, n_vars * n_lags, ifelse(check_roots, 1000, 1)))
  while(stationarity_check == FALSE) {
    iter <- iter + 1
    Pi_temp[,,iter] <- rmatn(M = t(post_pi), Q = post_pi_omega, P = Sigma_r)
    Pi_comp    <- build_companion(Pi_temp[,, iter], n_vars = n_vars, n_lags = n_lags)
    if (check_roots == TRUE) {
      root <- max_eig_cpp(Pi_comp)
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
#' @inherit pi_sigma_posterior
#' @templateVar Pi_r TRUE
#' @templateVar Sigma_r TRUE
#' @templateVar prior_psi TRUE
#' @templateVar prior_psi_omega TRUE
#' @templateVar D_mat TRUE
#' @templateVar n_determ TRUE
#' @template man_template
#' @return \code{psi_posterior} returns:
#' \item{psi_r}{The draw of \code{psi}.}
psi_posterior <- function(Pi_r, Sigma_r, Z_r1, prior_psi, prior_psi_omega, D_mat, n_vars, n_lags, n_determ) {
  U <- build_U_cpp(Pi = Pi_r, n_determ = n_determ,
                   n_vars = n_vars, n_lags = n_lags)
  post_psi_omega <- posterior_psi_omega(U = U, D_mat = D_mat, Sigma = Sigma_r,
                                        prior_psi_omega = prior_psi_omega)
  Y_tilde <- build_Y_tilde(Pi = Pi_r, z = Z_r1)

  post_psi <- posterior_psi(U = U, D_mat = D_mat, Sigma = Sigma_r, prior_psi_omega = prior_psi_omega,
                            post_psi_omega = post_psi_omega, Y_tilde = Y_tilde, prior_psi = prior_psi)
  psi_r <- t(rmultn(m = post_psi, Sigma = post_psi_omega))
  return(psi_r)
}

#' Draw from posterior of Z
#'
#' Function for drawing from the posterior of Z using a simulation smoother, which can be used as a block in a Gibbs sampler.
#' @inherit psi_posterior
#' @inherit pi_sigma_posterior
#' @templateVar Y TRUE
#' @templateVar psi_r TRUE
#' @templateVar Z_1 TRUE
#' @templateVar Lambda TRUE
#' @templateVar n_T_ TRUE
#' @templateVar smooth_state TRUE
#' @template man_template
#' @return \code{Z_posterior} returns a list with:
#' \item{Z_r}{The draw of \code{Z}.}
#' \item{smoothed_Z_r}{(Only if \code{smooth_state == TRUE}) The smoothed state.}
Z_posterior <- function(Y, d, Pi_r, Sigma_r, psi_r, Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state) {
  Q_comp     <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(Sigma_r))

  # Demean before putting into simulation smoother using the most recent draw of psi
  mZ <- Y - d %*% t(matrix(psi_r, nrow = n_vars))
  mZ <- as.matrix(mZ[-(1:n_lags), ])
  demeaned_z0 <- Z_1 - d[1:n_lags, ] %*% t(matrix(psi_r, nrow = n_vars))
  h0 <- matrix(t(demeaned_z0[n_lags:1,]), ncol = 1)
  Pi_comp <- build_companion(Pi_r, n_vars = n_vars, n_lags = n_lags)

  simulated_Z <- simulation_smoother(mZ = mZ, Lambda = Lambda, mF = Pi_comp, mQ = Q_comp, iT = n_T_,
                                     ip = n_vars, iq  = n_lags * n_vars, h0 = h0, P0 = diag(0, n_lags * n_vars))
  # For now, I'm just inserting h0 in the beginning of Z. Right now, Z has n_T number of rows,
  # but smooth_samp puts out n_T - n_lags rows.
  Z_r <- rbind(demeaned_z0, simulated_Z[,1:n_vars]) +
    d %*% t(matrix(psi_r, nrow = n_vars))

  # Also save the smoothed value of the state
  if (smooth_state == TRUE) {
    smoothed_state <- smoother(mZ = mZ, Lambda = Lambda, mF = Pi_comp, mQ = Q_comp, iT = n_T_,
                               ip = n_vars, iq  = n_lags * n_vars, h0 = h0, P0 = diag(0, n_lags * n_vars))
    smoothed_Z_r <- rbind(demeaned_z0, smoothed_state[, 1:n_vars]) +
      d %*% t(matrix(psi_r, nrow = n_vars))
    ret <- list(Z_r = Z_r, smoothed_Z_r = smoothed_Z_r)
  } else {
    ret <- list(Z_r = Z_r)
  }
  return(ret)
}

#' Compute posterior moments of the steady-state parameters
#'
#' Computes the mean and variance of the conditional posterior distribution of the steady-state parameters.
#' @templateVar U TRUE
#' @templateVar D_mat TRUE
#' @templateVar Sigma TRUE
#' @templateVar prior_psi_omega TRUE
#' @templateVar post_psi_omega TRUE
#' @templateVar Y_tilde TRUE
#' @templateVar prior_psi TRUE
#' @template man_template
#' @return The return is:
#' \item{psi}{The posterior mean (from \code{\link{posterior_psi}})}
posterior_psi <- function(U, D_mat, Sigma, prior_psi_omega, post_psi_omega, Y_tilde, prior_psi) {
  sigmaYD <- matrix(c(solve(Sigma) %*% t(Y_tilde) %*% D_mat), ncol = 1)
  psi <- post_psi_omega %*% (t(U) %*% sigmaYD + solve(prior_psi_omega) %*% prior_psi)
  return(psi)
}

#' @rdname posterior_psi
#' @return \item{psi_omega}{The posterior variance (from \code{\link{posterior_psi_omega}})}
posterior_psi_omega <- function(U, D_mat, Sigma, prior_psi_omega) {
  psi_omega <- solve(t(U) %*% (kronecker(crossprod(D_mat), solve(Sigma))) %*% U + solve(prior_psi_omega))
  return(psi_omega)
}


eval_Pi_Sigma_RaoBlack <- function(Z, d, post_psi, post_pi_mean, post_Sigma, nu, prior_pi, prior_pi_omega, prior_s, n_vars, n_lags, n_reps) {
  ################################################################
  ### Compute the Rao-Blackwellized estimate of posterior

  evals <- vector("numeric", n_reps - 1)

  for (i in 1:length(evals)) {
    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,,i+1] - d %*% post_psi
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    pi_sample <- solve(crossprod(XX)) %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_pi_omega_i <- solve(solve(prior_pi_omega) + crossprod(XX))
    post_pi_i       <- post_pi_omega_i %*% (solve(prior_pi_omega) %*% prior_pi + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% pi_sample)
    pi_diff <- prior_pi - pi_sample
    post_s_i <- prior_s + s_sample + t(pi_diff) %*% solve(post_pi_omega_i + solve(crossprod(XX))) %*% pi_diff

    # Evaluate
    evals[i] <- dnorminvwish(X = t(post_pi_mean), Sigma = post_Sigma, M = post_pi_i, P = post_pi_omega_i, S = post_s_i, v = nu)
  }

  return(evals)

}

eval_psi_MargPost <- function(Pi, Sigma, Z, post_psi, prior_psi, prior_psi_omega, D, n_determ, n_vars, n_lags, n_reps) {
  post_psi <- matrix(post_psi, ncol = 1)

  evals <- vector("numeric", n_reps - 1)

  ################################################################
  ### Steady-state step
  for (r in 1:(n_reps - 1)) {
    U <- build_U_cpp(Pi = Pi[,,r], n_determ = n_determ,
                     n_vars = n_vars, n_lags = n_lags)
    post_psi_omega <- posterior_psi_omega(U = U, D_mat = D, Sigma = Sigma[,, r],
                                          prior_psi_omega = prior_psi_omega)
    Y_tilde <- build_Y_tilde(Pi = Pi[,, r], z = Z[,, r])

    post_psi <- posterior_psi(U = U, D_mat = D, Sigma = Sigma[,, r], prior_psi_omega = prior_psi_omega,
                              post_psi_omega = post_psi_omega, Y_tilde = Y_tilde, prior_psi = prior_psi)

    evals[r] <- dmultn(x = post_psi, m = post_psi, Sigma = post_psi_omega)
  }

  return(evals)
}
