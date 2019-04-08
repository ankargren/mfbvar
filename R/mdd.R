#' Marginal data density estimation
#'
#' \code{mdd} estimates the (log) marginal data density.
#'
#' This is a generic function. See the methods for more information.
#' @seealso \code{\link{mdd.mfbvar_ss_iw}}, \code{\link{mdd.mfbvar_minn_iw}}
#' @param x argument to dispatch on (of class \code{mfbvar_ss} or \code{mfbvar_minn})
#' @param ... additional named arguments passed on to the methods

mdd <- function(x, ...) {
  UseMethod("mdd")
}

mdd.default <- function(x, ...) {
  stop("The marginal data density can currently only be estimated when inverse Wishart is used for the error covariance matrix.")
}

#' Marginal data density method for class \code{mfbvar_ss}
#'
#' Estimate the marginal data density for the model with a steady-state prior.
#' @param x object of class \code{mfbvar_ss}
#' @param method option for which method to choose for computing the mdd (\code{1} or \code{2})
#' @param ... additional arguments (currently only \code{p_trunc} for the degree of truncation for method 2 is available)
#' @details Two methods for estimating the marginal data density are implemented. Method 1 and 2 correspond to the two methods proposed by
#' Fuentes-Albero and Melosi (2013) and Ankargren, Unosson and Yang (2018).
#' @return The logarithm of the marginal data density.
#' @references Fuentes-Albero, C. and Melosi, L. (2013) Methods for Computing Marginal Data Densities from the Gibbs Output.
#' \emph{Journal of Econometrics}, 175(2), 132-141, \url{https://doi.org/10.1016/j.jeconom.2013.03.002}\cr
#'  Ankargren, S., Unosson, M., & Yang, Y. (2018) A Mixed-Frequency Bayesian Vector Autoregression with a Steady-State Prior. Working Paper, Department of Statistics, Uppsala University No. 2018:3.
#' @seealso \code{\link{mdd}}, \code{\link{mdd.mfbvar_minn_iw}}
mdd.mfbvar_ss_iw <- function(x, method = 1, ...) {
  if (method == 1) {
    mdd_est <- estimate_mdd_ss_1(x)
  } else if (method == 2) {
    mdd_est <- estimate_mdd_ss_2(x, ...)
  } else {
    stop("method: Must be 1 or 2.")
  }
  return(c(mdd_est$log_mdd))
}

#' Marginal data density method for class \code{mfbvar_minn}
#'
#' Estimate the marginal data density for the model with a Minnesota prior.
#' @param x object of class \code{mfbvar_minn}
#' @param ... additional arguments (currently only \code{p_trunc} for the degree of truncation is available)
#' @return The logarithm of the marginal data density.
#' @details The method used for estimating the marginal data density is the proposal made by
#' Schorfheide and Song (2015).
#' @references
#' Schorfheide, F., & Song, D. (2015) Real-Time Forecasting With a Mixed-Frequency VAR. \emph{Journal of Business & Economic Statistics}, 33(3), 366--380. \url{http://dx.doi.org/10.1080/07350015.2014.954707}
#' @seealso \code{\link{mdd}}, \code{\link{mdd.mfbvar_ss_iw}}
mdd.mfbvar_minn_iw <- function(x, ...) {
  quarterly_cols <- which(x$mfbvar_prior$freq == "q")
  estimate_mdd_minn(x, ...)
}

#' Estimate marginal data density in steady-state MF-BVAR
#'
#' This function provides the possibility to estimate the log marginal density using the steady-state MF-BVAR.
#' @keywords internal
#' @return
#' \code{estimate_mdd_ss_1} returns a list with components (all are currently in logarithms):
#' \item{lklhd}{The likelihood.}
#' \item{eval_prior_Pi_Sigma}{The evaluated prior.}
#' \item{eval_prior_psi}{The evaluated prior of psi.}
#' \item{eval_RB_Pi_Sigma}{The Rao-Blackwellized estimate of the conditional posterior of Pi and Sigma.}
#' \item{eval_marg_psi}{The evaluated marginal posterior of psi.}
#' \item{log_mdd}{The mdd estimate (in log).}
estimate_mdd_ss_1 <- function(mfbvar_obj) {
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
  Z     <- mfbvar_obj$Z
  d     <- mfbvar_obj$d
  Pi    <- mfbvar_obj$Pi
  Sigma <- mfbvar_obj$Sigma

  Lambda <- mfbvar_obj$Lambda

  post_Pi_mean <- apply(Pi, c(1, 2), mean)
  post_Sigma <- apply(Sigma, c(1, 2), mean)
  post_psi <- colMeans(psi)

  prior_S <- mfbvar_obj$prior_S
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_psi_Omega <- mfbvar_obj$prior_psi_Omega
  prior_psi_mean <- mfbvar_obj$prior_psi_mean

  freq <- mfbvar_obj$mfbvar_prior$freq
  Lambda <- build_Lambda(freq, n_lags)
  n_q <- sum(freq == "q")
  T_b <- max(which(!apply(apply(Y[, freq == "m"], 2, is.na), 1, any)))
  Lambda_ <- build_Lambda(rep("q", n_q), 3)

  ################################################################
  ### Initialize
  Pi_red    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma_red <- array(NA, dim = c(n_vars, n_vars, n_reps))
  Z_red <- array(NA, dim = c(n_T, n_vars, n_reps))

  Pi_red[,, 1]    <- post_Pi_mean
  Sigma_red[,, 1] <- post_Sigma
  Z_red[,, 1] <- apply(mfbvar_obj$Z, c(1, 2), mean)

  roots <- vector("numeric", n_reps)
  num_tries <- roots

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  D <- build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_Pi_Omega <- solve(prior_Pi_Omega)
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  Z_1 <- Z_red[1:n_lags,, 1]

  mZ <- Y - d %*% t(matrix(post_psi, nrow = n_vars))
  mZ <- as.matrix(mZ)
  demeaned_z0 <- Z_1 - d[1:n_lags, ] %*% t(matrix(post_psi, nrow = n_vars))
  d_post_psi <- d %*% t(matrix(post_psi, nrow = n_vars))
  ################################################################
  ### Reduced Gibbs step
  for (r in 2:n_reps) {
    ################################################################
    ### Pi and Sigma step
    #                             (Z_r1,                 d,     psi_r1,            prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi_Sigma <- posterior_Pi_Sigma(Z_r1 = Z_red[,, r-1], d = d, psi_r1 = post_psi, prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, n_vars+2, check_roots = TRUE, n_vars, n_lags, n_T)
    Pi_red[,,r]      <- Pi_Sigma$Pi_r
    Sigma_red[,,r]   <- Pi_Sigma$Sigma_r
    num_tries[r] <- Pi_Sigma$num_try
    roots[r]     <- Pi_Sigma$root

    ################################################################
    ### Smoothing step
    #(Y, d, Pi_r,            Sigma_r,               psi_r,                          Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state)

    Pi_r <- cbind(Pi_red[,,r], 0)
    Z_res <- kf_sim_smooth(mZ, Pi_r, Sigma_red[,,r], Lambda_, demeaned_z0, n_q, T_b)
    Z_res <- rbind(demeaned_z0, Z_res) + d_post_psi
    Z_red[,, r] <- Z_res
  }

  ################################################################
  ### For the likelihood calculation
  mZ <- Y - d %*% t(matrix(post_psi, nrow = n_vars))
  mZ <- mZ[-(1:n_lags), ]
  demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(post_psi, nrow = n_vars))
  h0 <- matrix(t(demeaned_z0), ncol = 1)
  h0 <- h0[(n_vars*n_lags):1,, drop = FALSE] # have to reverse the order
  Pi_comp <- build_companion(post_Pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp  <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0 <- matrix(0, n_lags*n_vars, n_lags*n_vars)

  ################################################################
  ### Final calculations
  lklhd          <- sum(c(loglike(Y = as.matrix(mZ), Lambda = Lambda, Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_, n_vars = n_vars, n_comp = n_lags * n_vars, z0 = h0, P0 = P0)[-1]))
  eval_prior_Pi_Sigma <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = prior_Pi_mean, P = prior_Pi_Omega, S = prior_S, v = n_vars+2)
  eval_prior_psi      <- dmultn(x = post_psi, m = prior_psi_mean, Sigma = prior_psi_Omega)
  eval_RB_Pi_Sigma    <- log(mean(eval_Pi_Sigma_RaoBlack(Z_array = Z_red, d = d, post_psi_center = post_psi, post_Pi_center = post_Pi_mean, post_Sigma_center = post_Sigma,
                                                         post_nu = post_nu, prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S,
                                                         n_vars = n_vars, n_lags = n_lags, n_reps = n_reps)))
  eval_marg_psi   <- log(mean(eval_psi_MargPost(Pi_array = Pi, Sigma_array = Sigma, Z_array = Z, post_psi_center = post_psi, prior_psi_mean = prior_psi_mean,
                                                prior_psi_Omega = prior_psi_Omega, D_mat = D, n_determ = n_determ, n_vars = n_vars, n_lags = n_lags, n_reps = n_reps)))

  mdd_estimate <- c(lklhd + eval_prior_Pi_Sigma + eval_prior_psi - (eval_RB_Pi_Sigma + eval_marg_psi))

  return(list(lklhd = lklhd, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma, eval_prior_psi = eval_prior_psi, eval_RB_Pi_Sigma = eval_RB_Pi_Sigma, eval_marg_psi = eval_marg_psi, log_mdd = mdd_estimate))
}


#' @rdname estimate_mdd_ss_1
#' @details \code{estimate_mdd_ss_1} uses method 1, \code{estimate_mdd_ss_2} uses method 2.
#' @templateVar mfbvar_obj TRUE
#' @templateVar p_trunc TRUE
#' @template man_template
#' @keywords internal
#' @return
#' \code{estimate_mdd_ss_1} returns a list with components being \code{n_reps}-long vectors and a scalar (the final estimate).
#' \item{eval_posterior_Pi_Sigma}{Posterior of Pi and Sigma.}
#' \item{data_likelihood}{The likelihood.}
#' \item{eval_prior_Pi_Sigma}{Prior of Pi and Sigma.}
#' \item{eval_prior_psi}{Prior of psi.}
#' \item{psi_truncated}{The truncated psi pdf.}
#' \item{log_mdd}{The mdd estimate (in log).}

estimate_mdd_ss_2 <- function(mfbvar_obj, p_trunc) {
  # Get things from the MFBVAR object
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

  Y <- mfbvar_obj$Y
  Z <- mfbvar_obj$Z
  d <- mfbvar_obj$d

  Lambda <- mfbvar_obj$Lambda

  post_Pi_mean <- apply(mfbvar_obj$Pi, c(1, 2), mean)
  post_Sigma <- apply(mfbvar_obj$Sigma, c(1, 2), mean)
  post_psi <- colMeans(psi)
  post_psi_Omega <- cov(psi)

  prior_S <- mfbvar_obj$prior_S
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_psi_Omega <- mfbvar_obj$prior_psi_Omega
  prior_psi_mean <- mfbvar_obj$prior_psi_mean

  # For the truncated normal
  chisq_val <- qchisq(p_trunc, n_determ*n_vars)

  #(mZ,lH,mF,mQ,iT,ip,iq,h0,P0)
  Pi_comp <- build_companion(post_Pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp  <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0      <- matrix(0, n_lags*n_vars, n_lags*n_vars)

  eval_posterior_Pi_Sigma <- vector("numeric", n_reps)
  data_likelihood <- vector("numeric", n_reps)
  eval_prior_Pi_Sigma <- vector("numeric", 1)
  eval_prior_psi <- vector("numeric", n_reps)
  psi_truncated <- vector("numeric", n_reps)

  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean
  for (r in 1:n_reps) {
    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,, r] - d %*% t(matrix(psi[r, ], nrow = n_vars))
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

    # Set the variables which vary in the Kalman filtering
    mZ <- Y - d %*% t(matrix(psi[r, ], nrow = n_vars))
    mZ <- mZ[-(1:n_lags), ]
    demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(psi[r, ], nrow = n_vars))
    h0 <- matrix(t(demeaned_z0), ncol = 1)
    h0 <- h0[(n_vars*n_lags):1,,drop = FALSE] # have to reverse the order


    eval_posterior_Pi_Sigma[r] <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = post_Pi_i, P = post_Pi_Omega_i, S = post_s_i, v = post_nu)
    data_likelihood[r] <- sum(c(loglike(Y = as.matrix(mZ), Lambda = Lambda, Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_, n_vars = n_vars, n_comp = n_lags * n_vars, z0 = h0, P0 = P0)[-1]))

    eval_prior_psi[r] <- dmultn(x = psi[r, ], m = prior_psi_mean, Sigma = prior_psi_Omega)
    psi_truncated[r] <- dnorm_trunc(psi[r, ], post_psi, solve(post_psi_Omega), n_determ*n_vars, p_trunc, chisq_val)

  }

  eval_prior_Pi_Sigma <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = prior_Pi_mean, P = prior_Pi_Omega, S = prior_S, v = n_vars+2)
  exp_term <- eval_posterior_Pi_Sigma - (data_likelihood + eval_prior_Pi_Sigma + eval_prior_psi)
  log_mdd <- -mean(exp_term)-log(mean(exp(exp_term-mean(exp_term)) * psi_truncated))

  return(list(eval_posterior_Pi_Sigma = eval_posterior_Pi_Sigma, data_likelihood = data_likelihood, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma,
              eval_prior_psi = eval_prior_psi, psi_truncated = psi_truncated, log_mdd = log_mdd))
}


#' Estimate marginal data density in Minnesota MF-BVAR
#'
#' This function provides the possibility to estimate the log marginal density (up to a constant) using the Minnesota MF-BVAR.
#' @rdname mdd.minn
#' @templateVar mfbvar_obj TRUE
#' @template man_template
#' @param quarterly_cols numeric vector with positions of quarterly variables
#' @templateVar p_trunc TRUE
#' @keywords internal
#' @return The log marginal data density estimate (bar a constant)
#'
estimate_mdd_minn <- function(mfbvar_obj, p_trunc, ...) {
  Z <- mfbvar_obj$Z
  Y <- mfbvar_obj$Y
  n_T <- dim(Z)[1]
  n_reps <- dim(Z)[3]
  n_vars <- ncol(Y)
  n_lags <- mfbvar_obj$mfbvar_prior$n_lags
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  inv_prior_Pi_Omega <- chol2inv(chol(prior_Pi_Omega))
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean
  prior_S <- mfbvar_obj$prior_S
  prior_nu <- mfbvar_obj$prior_nu

  postsim <- sapply(1:n_reps, function(x) {
    Z_comp <- build_Z(z = Z[,, x], n_lags = n_lags)
    XX <- Z_comp[-nrow(Z_comp), ]
    XX <- cbind(XX, 1)
    YY <- Z_comp[-1, 1:n_vars]

    XXt.XX <- crossprod(XX)
    XXt.XX.inv <- chol2inv(chol(XXt.XX))
    Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)

    # Posterior moments of Pi
    post_Pi_Omega <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
    post_Pi       <- post_Pi_Omega %*% (Omega_Pi + crossprod(XX, YY))
    S <- crossprod(YY - XX %*% Pi_sample)
    Pi_diff <- prior_Pi_mean - Pi_sample
    post_S <- prior_S + S + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff
    return(dmatt(YY, XX %*% prior_Pi_mean, chol2inv(chol(diag(nrow(YY)) + XX %*% prior_Pi_Omega %*% t(XX))), prior_S, prior_nu))
  })

  temp <- apply(Z[-(1:n_lags), , ], 3, function(x) x[is.na(c(mfbvar_obj$Y[-(1:n_lags),]))])

  if (length(temp) == 0) {
    return(log_mdd = mean(postsim))
  } else {
    n_para <- nrow(temp)
    drawmean <- matrix(rowMeans(temp), ncol = 1)
    drawsig <- cov(t(temp))
    drawsiginv <- chol2inv(chol(drawsig))
    drawsiglndet <- as.numeric(determinant(drawsiginv, logarithm = TRUE)$modulus)
    paradev <- temp - kronecker(matrix(1, 1, n_reps), drawmean)
    quadpara <- rowSums((t(paradev) %*% drawsiginv) * t(paradev))
    pcrit <- qchisq(p_trunc, df = nrow(drawmean))
    invlike <- matrix(NA, n_reps, length(p_trunc))
    indpara <- invlike
    lnfpara <- indpara
    densfac <- -0.5 * n_para * log(2 * pi) + 0.5 * drawsiglndet -
      0.5 * quadpara[1] - log(p_trunc) - postsim[1]
    densfac <- -mean(densfac)
    for (i in seq_along(p_trunc)) {
      for (j in 1:n_reps) {
        lnfpara[j, i] <- -0.5 * n_para * log(2 * pi) +
          0.5 * drawsiglndet - 0.5 * quadpara[j] - log(p_trunc[i])
        indpara[j, i] <- quadpara[j] < pcrit[i]
        invlike[j, i] <- exp(lnfpara[j, i] - postsim[j] +
                               densfac) * indpara[j, i]
      }
      meaninvlike <- colMeans(invlike)
      mdd <- densfac - log(meaninvlike)
    }
    return(log_mdd = mean(mdd) + sum(!is.na(Y[-(1:n_lags), mfbvar_obj$freq == "q"]))*log(3))
  }
}

