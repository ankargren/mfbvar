mdd1_qf <- function(mfbvar_obj) {
  ################################################################
  ### Get things from the MFBVAR object
  n_determ <- mfbvar_obj$n_determ
  n_vars <- mfbvar_obj$n_vars
  n_lags <- mfbvar_obj$n_lags
  n_T <- mfbvar_obj$n_T
  n_T_ <- mfbvar_obj$n_T_
  n_reps <- mfbvar_obj$n_reps

  psi <- mfbvar_obj$psi
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_S <- mfbvar_obj$prior_S
  post_nu <- mfbvar_obj$post_nu

  Y     <- mfbvar_obj$Y
  d     <- mfbvar_obj$d
  Pi    <- mfbvar_obj$Pi
  Sigma <- mfbvar_obj$Sigma

  Lambda <- mfbvar_obj$Lambda

  post_Pi_mean <- apply(Pi, c(1, 2), mean)
  post_Sigma <- apply(Sigma, c(1, 2), mean)
  post_psi <- colMeans(psi)

  prior_S <- mfbvar_obj$prior_S
  prior_nu <- mfbvar_obj$prior_nu
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_psi_Omega <- mfbvar_obj$prior_psi_Omega
  prior_psi_mean <- mfbvar_obj$prior_psi_mean

  #(mZ,lH,mF,mQ,iT,ip,iq,h0,P0)
  Pi_comp <- build_companion(post_Pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp  <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0      <- matrix(0, n_lags*n_vars, n_lags*n_vars)

  ################################################################
  ### Initialize
  Pi_red    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma_red <- array(NA, dim = c(n_vars, n_vars, n_reps))
  Pi_red[,, 1]    <- post_Pi_mean
  Sigma_red[,, 1] <- post_Sigma

  roots <- vector("numeric", n_reps)
  num_tries <- roots

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  D <- build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_Pi_Omega <- solve(prior_Pi_Omega)
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  ################################################################
  ### For the likelihood calculation
  mZ <- Y - d %*% t(matrix(post_psi, nrow = n_vars))
  mZ <- mZ[-(1:n_lags), ]
  demeaned_z0 <- Y[1:n_lags,] - d[1:n_lags, ] %*% t(matrix(post_psi, nrow = n_vars))
  h0 <- matrix(t(demeaned_z0), ncol = 1)
  h0 <- h0[(n_vars*n_lags):1,, drop = FALSE] # have to reverse the order

  ################################################################
  ### Final calculations
  # Eq (17) in Fuentes-Albero and Melosi

  ####
  # Likelihood
  Lambda <- matrix(0, nrow = n_vars, ncol = n_vars*n_lags)
  lklhd          <- sum(c(loglike(Y = as.matrix(mZ), Lambda = Lambda, Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_, n_vars = n_vars, n_comp = n_lags * n_vars, z0 = h0, P0 = P0)[-1]))

  ####
  # p(Pi, Sigma)
  eval_prior_Pi_Sigma <- log(dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = prior_Pi_mean, P = prior_Pi_Omega, S = prior_S, v = prior_nu))

  ####
  # p(psi)
  eval_prior_psi      <- log(dmultn(x = post_psi, m = prior_psi_mean, Sigma = prior_psi_Omega))

  ####
  # p(Pi, Sigma|psi, Y)
  demeaned_z <- Y - d %*% post_psi
  demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
  XX <- demeaned_Z[-nrow(demeaned_Z), ]
  YY <- demeaned_Z[-1, 1:n_vars]
  XXt.XX <- crossprod(XX)
  XXt.XX.inv <- chol2inv(chol(XXt.XX))
  Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)

  # Posterior moments of Pi
  post_Pi_Omega_i <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
  post_Pi_i       <- post_Pi_Omega_i %*% (Omega_Pi + crossprod(XX, YY))

  # Then Sigma
  s_sample  <- crossprod(YY - XX %*% Pi_sample)
  Pi_diff <- prior_Pi_mean - Pi_sample
  post_s_i <- prior_S + s_sample + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff

  # Evaluate
  eval_RB_Pi_Sigma <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = post_Pi_i, P = post_Pi_Omega_i, S = post_s_i, v = post_nu)

  ####
  # p(psi|Y)
  Z_array <- array(Y, dim = c(dim(Y)[1], dim(Y)[2], dim(Pi)[3]))
  eval_marg_psi   <- log(mean(eval_psi_MargPost(Pi_array = Pi, Sigma_array = Sigma, Z_array = Z_array, post_psi_center = post_psi, prior_psi_mean = prior_psi_mean,
                                                prior_psi_Omega = prior_psi_Omega, D_mat = D, n_determ = n_determ, n_vars = n_vars, n_lags = n_lags, n_reps = n_reps)))

  mdd_estimate <- lklhd + eval_prior_Pi_Sigma + eval_prior_psi - (eval_RB_Pi_Sigma + eval_marg_psi)

  return(list(lklhd = lklhd, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma, eval_prior_psi = eval_prior_psi, eval_RB_Pi_Sigma = eval_RB_Pi_Sigma, eval_marg_psi = eval_marg_psi, log_mdd = mdd_estimate))
}
