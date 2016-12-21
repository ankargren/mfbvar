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
