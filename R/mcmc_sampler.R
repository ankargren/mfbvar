#' MCMC sampler
#'
#' \code{mcmc_sampler} is a generic function for deciding which specific MCMC
#' algorithm to dispatch to. See the methods for more information.
#' @seealso \code{\link{mcmc_sampler.mfbvar_ss}}, \code{\link{mcmc_sampler.mfbvar_minn}}
#' @param x argument to dispatch on (of class \code{prior_obj} or \code{prior_obj})
#' @param ... additional named arguments passed on to the methods

mcmc_sampler <- function(x, ...) {
  UseMethod("mcmc_sampler")
}

#' MCMC sampler for mixed-frequency BVAR
#'
#' MCMC sampler to approximate the posterior distribution of the VAR model
#' parameters.
#'
#' @details
#' \code{mcmc_sampler.mfbvar_ss} is used for a steady-state prior and
#' \code{mcmc_sampler.mfbvar_minn} for a Minnesota prior
#'
#' @param x a prior object (inhereting from \code{mfbvar_prior})
#' @param ... additional arguments (\code{n_reps} and \code{init})
#' @keywords internal
#'
#' @return
#' An object of class \code{mfbvar} and \code{mfbvar_ss} or \code{mfbvar_minn}.

mcmc_sampler.mfbvar_ss <- function(x, ...) {

  n_vars <- ncol(x$Y)
  if (!(!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$prior_psi_Omega) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements:", names(test_sub)[which(test_sub)])
  }
  if (x$n_fcst > 0 && nrow(x$d_fcst) != x$n_fcst) {
    stop("d_fcst has ", nrow(x$d_fcst), " rows, but n_fcst is ", x$n_fcst, ".")
  }

  priors <- prior_Pi_Sigma(lambda1 = x$lambda1, lambda2 = x$lambda2, prior_Pi_AR1 = x$prior_Pi_AR1, Y = x$Y,
                           n_lags = x$n_lags, prior_nu = n_vars + 2)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  Y <- x$Y
  d <- x$d
  d_fcst <- x$d_fcst
  freq <- x$freq
  prior_psi_mean <- x$prior_psi_mean
  prior_psi_Omega <- x$prior_psi_Omega
  n_fcst <- x$n_fcst
  smooth_state <- x$smooth_state
  check_roots <- x$check_roots
  verbose <- x$verbose

  add_args <- list(...)
  n_reps <- add_args$n_reps
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_psi <- init$init_psi
  init_Z <- init$init_Z

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)
  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_Pi_mean)))/n_vars^2
  Lambda <- build_Lambda(freq, n_lags)
  n_pseudolags <- dim(Lambda)[2]/n_vars
  n_determ <- dim(d)[2]
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags


  n_q <- sum(freq == "q")
  T_b <- max(which(!apply(apply(Y[, freq == "m", drop = FALSE], 2, is.na), 1, any)))
  Lambda_ <- build_Lambda(rep("q", n_q), 3)
  ################################################################
  ### Preallocation
  # Pi and Sigma store their i-th draws in the third dimension, psi
  # is vectorized so it has its i-th draw stored in the i-th row
  # Pi:    p * pk * n_reps, each [,,i] stores Pi'
  # Sigma: p * p  * n_reps
  # psi:   n_reps * p
  # Z:     T * p * n_reps
  ### If forecasting (h is horizon):
  # Z_fcst: hk * p * n_reps
  # d_fcst_lags: hk * m
  ### If root checking:
  # roots: n_reps vector
  # num_tries: n_reps vector
  ### If smoothing of the state vector:
  # smoothed_Z: T * p * n_reps

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps))
  psi   <- array(NA, dim = c(n_reps, n_vars * n_determ))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  if (n_fcst > 0) {
    Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                   dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
    d_fcst_lags <- matrix(rbind(d[(n_T-n_lags+1):n_T, , drop = FALSE], d_fcst), nrow = n_fcst + n_lags)
  }
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  }
  if (smooth_state == TRUE) {
    smoothed_Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  }



  ################################################################
  ### MCMC sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the MCMC sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  if (is.null(init_Z)) {
    Z[,, 1] <- fill_na(Y)
  } else {
    if (all(dim(Z[,, 1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,, 1]), collapse = " x ")))
    }

  }

  ols_results <- tryCatch(ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = n_determ),
                          error = function(cond) NULL)
  if (is.null(ols_results)) {
    ols_results <- list()
    ols_results$Pi <- prior_Pi_mean
    ols_results$S <- prior_S
    ols_results$psi <- prior_psi_mean
  }

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- ols_results$Pi
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,, 1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- max_eig_cpp(Pi_comp)
  }

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    if (all(dim(Sigma[,,1]) == dim(init_Sigma))) {
      Sigma[,, 1] <- init_Sigma
    } else {
      stop(paste0("The dimension of init_Sigma is ", paste(dim(init_Sigma), collapse = " x "), ", but should be ", paste(dim(Sigma[,,1]), collapse = " x ")))
    }
  }

  if (is.null(init_psi)) {
    if (roots[1] < 1) {
      psi[1, ] <- ols_results$psi
    } else {
      psi[1, ] <- prior_psi_mean
    }
  } else {
    if (length(psi[1, ]) == length(init_psi)) {
      psi[1,] <- init_psi
    } else {
      stop(paste0("The length of init_psi is ", paste(length(init_psi), collapse = " x "), ", but should be ", paste(length(psi[1,]), collapse = " x ")))
    }
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  # if requested
  D_mat <- build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  Z_1 <- Z[1:n_pseudolags,, 1]

  if (verbose == TRUE) {
    pb <- timerProgressBar(width = 35, char = "[=-]", style = 5)
  }

  for (r in 2:(n_reps)) {
    ################################################################
    ### Pi and Sigma step
    #(Z_r1,             d,     psi_r1,                            prior_Pi_mean, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi_Sigma <- posterior_Pi_Sigma(Z_r1 = Z[,, r-1], d = d, psi_r1 = psi[r-1, , drop = FALSE], prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, n_vars+2, check_roots, n_vars, n_lags, n_T)
    Pi[,,r]      <- Pi_Sigma$Pi_r
    Sigma[,,r]   <- Pi_Sigma$Sigma_r
    num_tries[r] <- Pi_Sigma$num_try
    roots[r]     <- Pi_Sigma$root

    ################################################################
    ### Steady-state step
    #(Pi_r,            Sigma_r,               Z_r1,             prior_psi_mean, prior_psi_Omega, D, n_vars, n_lags, n_determ)
    psi[r, ] <- posterior_psi(Pi_r = Pi[,, r], Sigma_r = Sigma[,, r], Z_r1 = Z[,, r-1], prior_psi_mean, prior_psi_Omega, D_mat, n_vars, n_lags, n_determ)

    ################################################################
    ### Smoothing step

    Pi_r <- cbind(Pi[,,r], 0)
    mZ <- Y - d %*% t(matrix(psi[r,], nrow = n_vars))
    mZ <- as.matrix(mZ)
    demeaned_z0 <- Z_1 - d[1:n_lags, ] %*% t(matrix(psi[r,], nrow = n_vars))
    Z_res <- kf_sim_smooth(mZ, Pi_r, Sigma[,,r], Lambda_, demeaned_z0, n_q, T_b)
    Z_res <- rbind(demeaned_z0, Z_res) + d %*% t(matrix(psi[r,], nrow = n_vars))
    if (smooth_state == TRUE) {
      Z_smooth <- kf_ragged(mZ, Pi_r, Sigma[,,r], Lambda_, n_q, T_b)$Z_tT[-c(1:n_lags), ]
      smoothed_Z[,, r] <- rbind(demeaned_z0, Z_smooth) + d %*% t(matrix(psi[r, ], nrow = n_vars))
    }
    Z[,, r] <- Z_res

    ################################################################
    ### Forecasting step
    if (n_fcst > 0) {

      # Forecast the process with mean subtracted
      Z_fcst[1:n_lags, , r] <- Z[(n_T - n_lags+1):n_T,, r] - d[(n_T - n_lags+1):n_T, ] %*% t(matrix(psi[r, ], nrow = n_vars))
      for (h in 1:n_fcst) {
        Z_fcst[n_lags + h, , r] <- Pi[,, r] %*% matrix(c(t(Z_fcst[(n_lags+h-1):h,, r])), ncol = 1) +
          rmultn(m = matrix(0, nrow = n_vars), Sigma = Sigma[,,r])
      }

      # Add the mean
      Z_fcst[, , r] <- Z_fcst[, , r] + d_fcst_lags %*% t(matrix(psi[r, ], nrow = n_vars))
    }

    if (verbose == TRUE) {
      setTimerProgressBar(pb, r/n_reps)
    }
  }
  if (verbose == TRUE) {
    close(pb)
  }


  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = psi, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, mdd = NULL, smoothed_Z = NULL, n_determ = n_determ,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = prior_Pi_Omega, prior_Pi_mean = prior_Pi_mean,
                     prior_S = prior_S, prior_nu = n_vars+2, post_nu = n_T + n_vars+2, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = prior_psi_Omega, prior_psi_mean = prior_psi_mean, n_reps = n_reps - 1, Lambda = Lambda,
                     init = list(init_Pi = Pi[,, n_reps], init_Sigma = Sigma[,, n_reps], init_psi = psi[n_reps, ], init_Z = Z[,, n_reps]))

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  if (smooth_state == TRUE) {
    return_obj$smoothed_Z <- smoothed_Z
  }

  return(return_obj)

}

#' @rdname mcmc_sampler.mfbvar_ss

mcmc_sampler.mfbvar_minn <- function(x, ...){

  if (!(!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps))) {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    stop("Missing elements:", names(test_sub)[which(test_sub)])
  }

  Y <- x$Y
  freq <- x$freq
  prior_Pi_AR1 <- x$prior_Pi_AR1
  n_fcst <- x$n_fcst
  smooth_state <- x$smooth_state
  check_roots <- x$check_roots
  verbose <- x$verbose
  n_lags <- x$n_lags
  lambda1 <- x$lambda1
  lambda2 <- x$lambda2
  lambda3 <- x$lambda3

  add_args <- list(...)
  n_reps <- add_args$n_reps
  init <- add_args$init
  init_Pi <- init$init_Pi
  init_Sigma <- init$init_Sigma
  init_Z <- init$init_Z

  # n_vars: number of variables
  # n_lags: number of lags
  # n_determ: number of deterministic variables
  # n_T: sample size (full sample)
  # n_T_: sample size (reduced sample)

  Lambda <- build_Lambda(freq, n_lags)
  n_q <- sum(freq == "q")
  T_b <- max(which(!apply(apply(Y[, freq == "m", drop = FALSE], 2, is.na), 1, any)))
  Lambda_ <- build_Lambda(rep("q", n_q), 3)

  n_vars <- dim(Y)[2]
  n_pseudolags <- dim(Lambda)[2]/n_vars
  n_T <- dim(Y)[1]# - n_lags
  n_T_ <- n_T - n_pseudolags
  d <- matrix(1, nrow = nrow(Y), ncol = 1)
  lnpY1 <- rep(NA, n_reps)

  ################################################################
  ### Preallocation
  # Pi and Sigma store their i-th draws in the third dimension, psi
  # is vectorized so it has its i-th draw stored in the i-th row
  # Pi:    p * pk * n_reps, each [,,i] stores Pi'
  # Sigma: p * p  * n_reps
  # psi:   n_reps * p
  # Z:     T * p * n_reps
  ### If forecasting (h is horizon):
  # Z_fcst: hk * p * n_reps
  # d_fcst_lags: hk * m
  ### If root checking:
  # roots: n_reps vector
  # num_tries: n_reps vector
  ### If smoothing of the state vector:
  # smoothed_Z: T * p * n_reps

  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags + 1, n_reps))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_reps))
  Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  if (n_fcst>0) {
    Z_fcst<- array(NA, dim = c(n_fcst+n_lags, n_vars, n_reps),
                   dimnames = list(c((n_T-n_lags+1):n_T, paste0("fcst_", 1:n_fcst)), NULL, NULL))
  }
  if (check_roots == TRUE) {
    roots <- vector("numeric", n_reps)
    num_tries <- roots
  }
  if (smooth_state == TRUE) {
    smoothed_Z     <- array(NA, dim = c(n_T, n_vars, n_reps))
  }



  ################################################################
  ### MCMC sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the MCMC sampler for a burn-in
  # period, then use the final draw of that as initialization
  # for multiple chains

  if (is.null(init_Z)) {
    Z[,, 1] <- fill_na(Y)
  } else {
    if (all(dim(Z[,, 1]) == dim(init_Z))) {
      Z[,, 1] <- init_Z
    } else {
      stop(paste0("The dimension of init_Z is ", paste(dim(init_Z), collapse = " x "), ", but should be ", paste(dim(Z[,, 1]), collapse = " x ")))
    }

  }

  ols_results <- ols_initialization(z = Z[,, 1], d = d, n_lags = n_lags, n_T = n_T, n_vars = n_vars, n_determ = 1)

  if (is.null(init_Pi)) {
    Pi[,, 1]    <- cbind(ols_results$Pi, ols_results$const)
  } else {
    if (all(dim(Pi[,, 1]) == dim(init_Pi))) {
      Pi[,, 1] <- init_Pi
    } else {
      stop(paste0("The dimension of init_Pi is ", paste(dim(init_Pi), collapse = " x "), ", but should be ", paste(dim(Pi[,, 1]), collapse = " x ")))
    }
  }

  # Compute the maximum eigenvalue of the initial Pi
  if (check_roots == TRUE) {
    Pi_comp    <- build_companion(Pi = Pi[,-ncol(Pi[,,1]), 1], n_vars = n_vars, n_lags = n_lags)
    roots[1]   <- max_eig_cpp(Pi_comp)
  }

  if (is.null(init_Sigma)) {
    Sigma[,, 1] <- ols_results$S
  } else {
    if (all(dim(Sigma[,,1]) == dim(init_Sigma))) {
      Sigma[,, 1] <- init_Sigma
    } else {
      stop(paste0("The dimension of init_Sigma is ", paste(dim(init_Sigma), collapse = " x "), ", but should be ", paste(dim(Sigma[,,1]), collapse = " x ")))
    }
  }

  ################################################################
  ### Compute terms which do not vary in the sampler

  Z_1 <- Z[1:n_pseudolags,, 1]

  ####################################################
  Y_bar <- colMeans(Y, na.rm = TRUE)
  s_bar <- sqrt(diag(prior_Pi_Sigma(0.2, 1, prior_Pi_AR1, Y, n_lags, n_vars + 2)$prior_S))

  dummy_size <- (n_lags + 1)*n_vars + 3
  breaks <- numeric(5)
  Y_dum <- matrix(0, nrow = dummy_size, ncol = n_vars)
  X_dum <- matrix(0, nrow = dummy_size, ncol = n_vars*n_lags + 1)
  n_XX <- ncol(X_dum)
  ## 1: AR(1) coefficients
  Y_dum[1:n_vars, ] <- diag(s_bar * prior_Pi_AR1)/lambda1
  breaks[1] <- n_vars

  ## 2: AR(2), ..., AR(p) coefficients
  X_dum[1:(n_vars*n_lags), 1:(n_vars*n_lags)] <- kronecker(diag((1:n_lags)^lambda2), diag(s_bar))/lambda1
  if (n_lags > 1) {
    breaks[2] <- breaks[1] + (n_lags - 1)*n_vars
  } else {
    breaks[2] <- breaks[1]
  }

  ## 3: Sigma
  Y_dum[(breaks[2]+1):(breaks[2]+n_vars), ] <- diag(s_bar)
  breaks[3] <- breaks[2] + n_vars

  ## 4: Intercept
  X_dum[breaks[3] + 1, n_XX] <- 1/lambda3

  #
  #   ## 3
  #   Y_dum[(breaks[2]+1):(breaks[2]+lambda3*n_vars), ] <- kronecker(matrix(1, nrow = lambda3, ncol = 1), pre_Sigma)
  #   breaks[3] <- breaks[2] + lambda3 * n_vars
  #
  #   ##
  #   lambda_mean <- lambda4 * Y_bar
  #   Y_dum[breaks[3] + 1, ] <- lambda_mean
  #   X_dum[breaks[3] + 1, ] <- cbind(kronecker(matrix(1, nrow = 1, ncol = n_lags), lambda_mean), lambda4)
  #   breaks[4] <- breaks[3] + 1
  #
  #   ##
  #   mu_mean <- diag(lambda5 * Y_bar)
  #   Y_dum[(breaks[4] + 1):(breaks[4] + n_vars), ] <- mu_mean
  #   X_dum[(breaks[5] + 1):(breaks[4] + n_vars), 1:(n_lags * n_vars)] <- kronecker(matrix(1, nrow = 1, ncol = n_lags), mu_mean)
  #   breaks[5] <- breaks[4] + n_vars


  Pi_r1 <- Pi[,,1]
  const_r1 <- Pi_r1[, ncol(Pi_r1)]
  Pi_r1 <- Pi_r1[, -ncol(Pi_r1)]


  ####################################################

  # MDD
  # svd_res <- svd(crossprod(X_dum), nu = n_XX)
  # ux <- svd_res$u
  # sx <- rbind(diag(svd_res$d), matrix(0, nrow = dim(ux)[1]- length(svd_res$d), ncol = length(svd_res$d)))
  # vx <- svd_res$v
  # sv_XX <- ux %*% sqrt(sx) %*% t(vx)
  #
  # upx <- ux
  # spx <- sx
  # vpx <- vx
  # inv_spx <- matrix(0, ncol = nrow(spx), nrow = nrow(spx))
  #
  # for (rr in 1:nrow(spx)) {
  #   if (spx[rr, rr] > 1e-12) {
  #     inv_spx[rr, rr] <- 1/spx[rr, rr]
  #   }
  # }
  #
  # inv_sv_XX <- t(upx %*% inv_spx %*% t(vpx))


  n_dummy <- nrow(X_dum)
  Phi <- tcrossprod(chol2inv(chol(crossprod(X_dum))), crossprod(Y_dum, X_dum))
  S0 <- crossprod(Y_dum - X_dum %*% Phi)
  gam0 <- sum(lgamma(0.5 * (n_dummy - n_XX + 1 - 1:n_vars)))

  lnpY0 <- - n_vars * 0.5 * determinant(crossprod(X_dum), logarithm = TRUE)$modulus -
    (n_dummy - n_XX)*0.5*determinant(S0, logarithm = TRUE)$modulus + n_vars * (n_vars - 1) * 0.25*log(pi) + gam0


  if (verbose == TRUE) {
    pb <- timerProgressBar(width = 35, char = "[=-]", style = 5)
  }

  for (r in 2:(n_reps)) {

    Pi_r1 <- Pi[,,r-1]
    Sigma_r1 <- Sigma[,,r-1]

    Z_res <- kf_sim_smooth(Y, Pi_r1, Sigma_r1, Lambda_, Z_1, n_q, T_b)

    Z[,, r] <- rbind(Z_1, Z_res)

    Z_comp <- build_Z(z = Z[,, r], n_lags = n_lags)
    XX_act <- Z_comp[-nrow(Z_comp), ]
    XX_act <- cbind(XX_act, 1)
    YY_act <- Z_comp[-1, 1:n_vars]
    YY <- rbind(Y_dum, YY_act)
    XX <- rbind(X_dum, XX_act)

    XXt.XX <- crossprod(XX)
    XXt.XX.inv <- chol2inv(chol(XXt.XX))
    Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)
    post_Pi <- Pi_sample
    post_Pi_Omega <- XXt.XX.inv
    S <- crossprod(YY - XX %*% Pi_sample)

    post_nu <- nrow(YY) - ncol(YY)*n_vars - 1

    Sigma_r <- rinvwish(v = post_nu, S = S)
    Sigma[,, r]   <- Sigma_r

    Pi[,, r] <- rmatn(M = t(post_Pi), Q = post_Pi_Omega, P = Sigma_r)

    if (check_roots) {
      Pi_comp  <- build_companion(Pi[,-ncol(Pi_r1), r], n_vars = n_vars, n_lags = n_lags)
      roots[r] <- max_eig_cpp(Pi_comp)
    }


    # # MDD
    # svd_res <- svd(crossprod(XX), nu = ncol(XX))
    # ux <- svd_res$u
    # sx <- rbind(diag(svd_res$d), matrix(0, nrow = dim(ux)[1]- length(svd_res$d), ncol = length(svd_res$d)))
    # vx <- svd_res$v
    # sv_XX <- ux %*% sqrt(sx) %*% t(vx)
    #
    # upx <- ux
    # spx <- sx
    # vpx <- vx
    # inv_spx <- matrix(0, ncol = nrow(spx), nrow = nrow(spx))
    #
    # for (rr in 1:nrow(spx)) {
    #   if (spx[rr, rr] > 1e-12) {
    #     inv_spx[rr, rr] <- 1/spx[rr, rr]
    #   }
    # }
    #
    # inv_sv_XX <- t(upx %*% inv_spx %*% t(vpx))



    n_tot <- nrow(XX)
    Phi <- tcrossprod(chol2inv(chol(crossprod(XX))), crossprod(YY, XX))
    S1 <- crossprod(YY - XX %*% Phi)
    gam1 <- sum(lgamma(0.5 * (n_tot - n_XX + 1 - 1:n_vars)))

    lnpY1[r] <- - n_vars * 0.5 * determinant(XXt.XX, logarithm = TRUE)$modulus -
      (n_tot - n_XX)*0.5*determinant(S1, logarithm = TRUE)$modulus + n_vars * (n_vars - 1) * 0.25*log(pi) + gam1

    Pi_r <- Pi[,,r]
    const_r <- Pi_r[, ncol(Pi_r)]
    Pi_r <- Pi_r[, -ncol(Pi_r)]


    ################################################################
    ### Forecasting step
    if (n_fcst>0) {

      # Forecast the process with mean subtracted
      Z_fcst[1:n_lags, , r] <- Z[(n_T - n_lags+1):n_T,, r]
      for (h in 1:n_fcst) {
        Z_fcst[n_lags + h, , r] <- const_r + Pi_r  %*% matrix(c(t(Z_fcst[(n_lags+h-1):h,, r])), ncol = 1) +
          rmultn(m = matrix(0, nrow = n_vars), Sigma = Sigma[,,r])
      }

    }
    #########################################
    # Add the likelihood here
    if (verbose == TRUE) {
      setTimerProgressBar(pb, r/n_reps)
    }
  }

  if (verbose == TRUE) {
    close(pb)
  }


  ################################################################
  ### Prepare the return object
  return_obj <- list(Pi = Pi, Sigma = Sigma, psi = NULL, Z = Z, roots = NULL, num_tries = NULL,
                     Z_fcst = NULL, mdd = NULL, smoothed_Z = NULL, n_determ = 1,
                     n_lags = n_lags, n_vars = n_vars, n_fcst = n_fcst, prior_Pi_Omega = NULL, prior_Pi_mean = NULL,
                     prior_S = NULL, prior_nu = NULL, post_nu = NULL, d = d, Y = Y, n_T = n_T, n_T_ = n_T_,
                     prior_psi_Omega = NULL, prior_psi_mean = NULL, n_reps = n_reps-1, Lambda = Lambda,
                     lnpYY = lnpY1 - lnpY0, freq = freq,
                     init = list(init_Pi = Pi[,, n_reps], init_Sigma = Sigma[,, n_reps], init_Z = Z[,, n_reps]))

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }
  if (n_fcst>0) {
    return_obj$Z_fcst <- Z_fcst
  }
  if (smooth_state == TRUE) {
    return_obj$smoothed_Z <- smoothed_Z
  }

  return(return_obj)

}


