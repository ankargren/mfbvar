fill_na <- function(Y) {
  apply(Y, 2, function(x) {
    n_x <- length(x) # save lentgh
    if (any(is.na(x))) {
      x <- x[1:max(which(is.na(x) == FALSE))] # get rid of NAs in the end
      for (i in which(is.na(x))) {
        x1 <- NA
        counter <- 1
        while (is.na(x1) == TRUE) {
          x1 <- x[i + counter]
          counter <- counter + 1
        }
        x[i] <- x1
      }

      trimmed_length <- length(x)
      if (trimmed_length < n_x) {
        x <- c(x, rep(NA, n_x - trimmed_length))
        for (i in trimmed_length:n_x) {
          x[i] <- x[trimmed_length]
        }
      }
    }
    x})
}

dnorminvwish <- function(X, Sigma, M, P, S, v) {
  q <- dim(Sigma)[1]
  p <- dim(P)[1]
  det_Sigma <- det(Sigma)
  inv_Sigma <- solve(Sigma)
  dmultnorm <- (-p*q/2) * log(2 * pi) + (-p/2) * log(det_Sigma) + (-q/2)*log(det(P)) + (-1/2 * sum(diag(inv_Sigma %*% t(X - M) %*% solve(P) %*% (X - M))))
  cc <- (v * q/2)*log(2) + (q*(q-1)/4)*log(pi) + sum(lgamma((v+1-1:q)/2))
  dinvwish <- -cc + (v/2) * log(det(S)) -(v+q+1)/2*log(det_Sigma) -1/2 * sum(diag(inv_Sigma %*% S))
  return(exp(dmultnorm + dinvwish))
}

dmultn <- function(x, m, Sigma) {
  p <- dim(Sigma)[1]
  log_d <- (-1/2)* log(det(2*pi*Sigma)) -1/2 * t(x-m) %*% solve(Sigma) %*% (x-m)
  return(exp(log_d))
}

dnorm_trunc <- function(x, m, V_inv, d, p_trunc, chisq_val) {
  qf <- t(x - m) %*% V_inv %*% (x - m)
  return((1/p_trunc) * (1/sqrt((2*pi)^d/det(V_inv))) * exp(-0.5 * qf) * (qf < chisq_val))
}

mdd <- function(mfbvar_obj, p_trunc) {
  # Get things from the MFBVAR object
  n_determ <- mfbvar_obj$n_determ
  n_vars <- mfbvar_obj$n_vars
  n_lags <- mfbvar_obj$n_lags
  n_T <- mfbvar_obj$n_T
  n_T_ <- mfbvar_obj$n_T_
  n_reps <- mfbvar_obj$n_reps

  psi <- mfbvar_obj$psi
  prior_pi_omega <- mfbvar_obj$prior_pi_omega
  prior_pi <- mfbvar_obj$prior_pi
  prior_s <- mfbvar_obj$prior_s
  nu <- mfbvar_obj$nu

  Y <- mfbvar_obj$Y
  Z <- mfbvar_obj$Z
  d <- mfbvar_obj$d

  Lambda <- mfbvar_obj$Lambda

  post_pi_mean <- apply(mfbvar_obj$Pi, c(1, 2), mean)
  post_Sigma <- apply(mfbvar_obj$Sigma, c(1, 2), mean)
  post_psi <- colMeans(psi)
  post_psi_omega <- cov(psi)

  prior_s <- mfbvar_obj$prior_s
  prior_nu <- mfbvar_obj$prior_nu
  prior_pi_omega <- mfbvar_obj$prior_pi_omega
  prior_pi <- mfbvar_obj$prior_pi
  prior_psi_omega <- mfbvar_obj$prior_psi_omega
  prior_psi <- mfbvar_obj$prior_psi

  # For the truncated normal
  chisq_val <- qchisq(p_trunc, n_determ*n_vars)

  #(mZ,lH,mF,mQ,iT,ip,iq,h0,P0)
  Pi_comp <- build_companion(post_pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp  <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0      <- matrix(0, n_lags*n_vars, n_lags*n_vars)

  pi_sigma_posterior <- vector("numeric", n_reps)
  data_likelihood <- vector("numeric", n_reps)
  pi_sigma_prior <- vector("numeric", n_reps)
  psi_prior <- vector("numeric", n_reps)
  psi_truncated <- vector("numeric", n_reps)
  for (r in 1:n_reps) {
    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,, r] - d %*% t(matrix(psi[r, ], nrow = n_vars))
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

    # Set the variables which vary in the Kalman filtering
    mZ <- Y - d %*% t(matrix(psi[r, ], nrow = n_vars))
    mZ <- mZ[-(1:n_lags), ]
    demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(psi[r, ], nrow = n_vars))
    h0 <- matrix(t(demeaned_z0), ncol = 1)
    h0 <- h0[(n_vars*n_lags):1,,drop = FALSE] # have to reverse the order


    pi_sigma_posterior[r] <- dnorminvwish(X = t(post_pi_mean), Sigma = post_Sigma, M = post_pi_i, P = post_pi_omega_i, S = post_s_i, v = nu)
    data_likelihood[r] <- exp(sum(c(loglike(mZ = as.matrix(mZ), Lambda = Lambda, mF = Pi_comp, mQ = Q_comp, iT = n_T_, ip = n_lags, iq = n_lags * n_vars, h0 = h0, P0 = P0)[-1])))
    pi_sigma_prior[r] <- dnorminvwish(X = t(post_pi_mean), Sigma = post_Sigma, M = prior_pi, P = prior_pi_omega, S = prior_s, v = prior_nu)
    psi_prior[r] <- dmultn(x = psi[r, ], m = prior_psi, Sigma = prior_psi_omega)
    psi_truncated[r] <- dnorm_trunc(psi[r, ], post_psi, solve(post_psi_omega), n_determ*n_vars, p_trunc, chisq_val)

  }
  return(list(pi_sigma_posterior = pi_sigma_posterior, data_likelihood = data_likelihood, pi_sigma_prior = pi_sigma_prior,
              psi_prior = psi_prior, psi_truncated = psi_truncated))
}

mdd1 <- function(mfbvar_obj) {
  ################################################################
  ### Get things from the MFBVAR object
  n_determ <- mfbvar_obj$n_determ
  n_vars <- mfbvar_obj$n_vars
  n_lags <- mfbvar_obj$n_lags
  n_T <- mfbvar_obj$n_T
  n_T_ <- mfbvar_obj$n_T_
  n_reps <- mfbvar_obj$n_reps

  psi <- mfbvar_obj$psi
  prior_pi_omega <- mfbvar_obj$prior_pi_omega
  prior_pi <- mfbvar_obj$prior_pi
  prior_s <- mfbvar_obj$prior_s
  nu <- mfbvar_obj$nu

  Y     <- mfbvar_obj$Y
  Z     <- mfbvar_obj$Z
  d     <- mfbvar_obj$d
  Pi    <- mfbvar_obj$Pi
  Sigma <- mfbvar_obj$Sigma

  Lambda <- mfbvar_obj$Lambda

  post_pi_mean <- apply(Pi, c(1, 2), mean)
  post_Sigma <- apply(Sigma, c(1, 2), mean)
  post_psi <- colMeans(psi)
  post_psi_omega <- cov(psi)

  prior_s <- mfbvar_obj$prior_s
  prior_nu <- mfbvar_obj$prior_nu
  prior_pi_omega <- mfbvar_obj$prior_pi_omega
  prior_pi <- mfbvar_obj$prior_pi
  prior_psi_omega <- mfbvar_obj$prior_psi_omega
  prior_psi <- mfbvar_obj$prior_psi

  #(mZ,lH,mF,mQ,iT,ip,iq,h0,P0)
  Pi_comp <- build_companion(post_pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp  <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0      <- matrix(0, n_lags*n_vars, n_lags*n_vars)

  ################################################################
  ### Initialize
  Pi_red    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma_red <- array(NA, dim = c(n_vars, n_vars, n_reps))
  Z_red <- array(NA, dim = c(n_T, n_vars, n_reps))

  Pi_red[,, 1]    <- post_pi_mean
  Sigma_red[,, 1] <- post_Sigma
  Z_red[,, 1] <- apply(mfbvar_obj$Z, c(1, 2), mean)

  roots <- vector("numeric", n_reps)
  num_tries <- roots

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  D <- build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_pi_omega <- solve(prior_pi_omega)
  omega_pi <- inv_prior_pi_omega %*% prior_pi

  Z_1 <- Z_red[1:n_lags,, 1]

  ################################################################
  ### Reduced Gibbs step
  for (r in 2:n_reps) {
    ################################################################
    ### Pi and Sigma step
    #(Z_r1,             d,     psi_r1,                            prior_pi, inv_prior_pi_omega, omega_pi, prior_s, prior_nu, check_roots, n_vars, n_lags, n_T)
    pi_sigma <- pi_sigma_posterior(Z_r1 = Z_red[,, r-1], d = d, psi_r1 = post_psi, prior_pi, inv_prior_pi_omega, omega_pi, prior_s, prior_nu, check_roots = TRUE, n_vars, n_lags, n_T)
    Pi_red[,,r]      <- pi_sigma$Pi_r
    Sigma_red[,,r]   <- pi_sigma$Sigma_r
    num_tries[r] <- pi_sigma$num_try
    roots[r]     <- pi_sigma$root

    ################################################################
    ### Smoothing step
    #(Y, d, Pi_r,            Sigma_r,               psi_r,                          Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state)
    Z_res <- Z_posterior(Y, d, Pi_r = Pi_red[,, r], Sigma_r = Sigma_red[,, r], psi_r = post_psi, Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state = FALSE)
    Z_red[,, r] <- Z_res$Z_r
  }

  ################################################################
  ### For the likelihood calculation
  mZ <- Y - d %*% t(matrix(post_psi, nrow = n_vars))
  mZ <- mZ[-(1:n_lags), ]
  demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(post_psi, nrow = n_vars))
  h0 <- matrix(t(demeaned_z0), ncol = 1)
  h0 <- h0[(n_vars*n_lags):1,, drop = FALSE] # have to reverse the order

  ################################################################
  ### Final calculations
  lklhd          <- exp(sum(c(loglike(mZ = as.matrix(mZ), Lambda = Lambda, mF = Pi_comp, mQ = Q_comp, iT = n_T_, ip = n_lags, iq = n_lags * n_vars, h0 = h0, P0 = P0)[-1])))
  Pi_Sigma_prior <- dnorminvwish(X = t(post_pi_mean), Sigma = post_Sigma, M = prior_pi, P = prior_pi_omega, S = prior_s, v = prior_nu)
  psi_prior      <- dmultn(x = post_psi, m = prior_psi, Sigma = prior_psi_omega)
  Pi_Sigma_RB    <- mean(eval_Pi_Sigma_RaoBlack(Z_red, d, post_pi_mean, post_Sigma, nu, post_psi, prior_pi, prior_pi_omega, prior_s, n_vars, n_lags, n_reps))
  psi_MargPost   <- mean(eval_psi_MargPost(Pi, Sigma, Z, post_psi, prior_psi_omega, D, n_determ, n_vars, n_lags, n_reps))

  mdd_estimate <- lklhd * Pi_Sigma_prior * psi_prior / (Pi_Sigma_RB * psi_MargPost)

  return(list(lklhd = lklhd, Pi_Sigma_prior = Pi_Sigma_prior, psi_prior = psi_prior, Pi_Sigma_RB = Pi_Sigma_RB, psi_MargPost = psi_MargPost, mdd_estimate = mdd_estimate))
}
