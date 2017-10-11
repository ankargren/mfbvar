#' Gibbs sampler for Mixed-Frequency BVAR
#'
#' \code{gibbs_sampler} runs a Gibbs sampler to approximate the posterior of the VAR model parameters.
#'
#' @templateVar Y TRUE
#' @templateVar freq TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar lambda3 TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_reps TRUE
#' @templateVar init_Pi TRUE
#' @templateVar init_Sigma TRUE
#' @templateVar init_Z TRUE
#' @templateVar smooth_state TRUE
#' @templateVar check_roots TRUE
#' @templateVar verbose TRUE
#' @template man_template
#' @keywords internal
#' @details
#' The prior covariance of \eqn{\Pi} given \eqn{\Sigma} is \eqn{\Sigma \otimes \Omega_\Pi}.
#'
#' The function \code{gibbs_sampler_qf} is a temporary version of \code{gibbs_sampler} which is customized for quarterly data (for improved speed). Its purpose is to accomomdate conditional forecasting. The function \code{gibbs_sampler2} uses an R implementation of the simulation smoother, whereas \code{gibbs_sampler} uses a C++ implementation.
#'
#' @return
#' An object of class mfbvar.

gibbs_sampler_minn <- function(Y, freq, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst = NULL, n_reps,
                                 init_Pi = NULL, init_Sigma = NULL, init_Z = NULL,
                                 smooth_state = FALSE, check_roots = TRUE, verbose = TRUE){

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
  ### Gibbs sampling initialization

  # If the initial values are not provided, the missing values in
  # Z are filled with the next observed value and Pi, Sigma and
  # psi are then computed using maximum likelihood

  # This allows the user to run the Gibbs sampler for a burn-in
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

  dummy_size <- 1 + (n_lags + 2)*n_vars
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

    Z_res <- kf_sim_smooth(Y, Pi_r1, Sigma_r1, Lambda_, Z_1, n_q, T_b)[-c(1:n_lags), ]

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
                     lnpYY = lnpY1 - lnpY0, freq = freq)

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

