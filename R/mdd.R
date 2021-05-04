#' Marginal data density estimation
#'
#' \code{mdd} estimates the (log) marginal data density.
#'
#' This is a generic function. See the methods for more information.
#' @seealso \code{\link{mdd.mfbvar_ss_iw}}, \code{\link{mdd.mfbvar_minn_iw}}
#' @param x argument to dispatch on (of class \code{mfbvar_ss} or \code{mfbvar_minn})
#' @param ... additional named arguments passed on to the methods
#' @return The logarithm of the marginal data density.
#' @details The marginal data density is also known as the marginal likelihood.

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
#' @param ... additional arguments (currently unused)
#' @details The marginal data density estimator follows Fuentes-Albero and Melosi (2013) and Ankargren, Unosson and Yang (2018).
#' @return The logarithm of the marginal data density.
#' @references Fuentes-Albero, C. and Melosi, L. (2013) Methods for Computing Marginal Data Densities from the Gibbs Output.
#' \emph{Journal of Econometrics}, 175(2), 132-141, \doi{10.1016/j.jeconom.2013.03.002}\cr
#'  Ankargren, S., Unosson, M., & Yang, Y. (2018) A Mixed-Frequency Bayesian Vector Autoregression with a Steady-State Prior. Working Paper, Department of Statistics, Uppsala University No. 2018:3.
#' @seealso \code{\link{mdd}}, \code{\link{mdd.mfbvar_minn_iw}}
mdd.mfbvar_ss_iw <- function(x, ...) {
  if (x$mfbvar_prior$aggregation != "average") {
    stop("The marginal data density can only be computed using intra-quarterly average aggregation.")
  }
  mdd_est <- estimate_mdd_ss(x)
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
#' Schorfheide, F., & Song, D. (2015) Real-Time Forecasting With a Mixed-Frequency VAR. \emph{Journal of Business & Economic Statistics}, 33(3), 366--380. \doi{10.1080/07350015.2014.954707}
#' @seealso \code{\link{mdd}}, \code{\link{mdd.mfbvar_ss_iw}}
mdd.mfbvar_minn_iw <- function(x, ...) {
  if (x$mfbvar_prior$aggregation != "average") {
    stop("The marginal data density can only be computed using intra-quarterly average aggregation.")
  }
  quarterly_cols <- which(x$mfbvar_prior$freq == "q")
  estimate_mdd_minn(x, ...)
}

#' Estimate marginal data density in steady-state MF-BVAR
#'
#' This function provides the possibility to estimate the log marginal density using the steady-state MF-BVAR.
#' @keywords internal
#' @noRd
#' @return
#' \code{estimate_mdd_ss} returns a list with components (all are currently in logarithms):
#' \item{lklhd}{The likelihood.}
#' \item{eval_prior_Pi_Sigma}{The evaluated prior.}
#' \item{eval_prior_psi}{The evaluated prior of psi.}
#' \item{eval_RB_Pi_Sigma}{The Rao-Blackwellized estimate of the conditional posterior of Pi and Sigma.}
#' \item{eval_marg_psi}{The evaluated marginal posterior of psi.}
#' \item{log_mdd}{The mdd estimate (in log).}
estimate_mdd_ss <- function(mfbvar_obj) {
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
  post_nu <- mfbvar_obj$prior_nu + mfbvar_obj$n_T_

  Y <- mfbvar_obj$Y
  Z <- mfbvar_obj$Z
  d <- mfbvar_obj$d
  Pi <- mfbvar_obj$Pi
  Sigma <- mfbvar_obj$Sigma

  Lambda <- mfbvar_obj$Lambda

  post_Pi_mean <- apply(Pi, c(1, 2), mean)
  post_Sigma <- apply(Sigma, c(1, 2), mean)
  post_psi <- apply(psi, c(1, 2), mean)

  prior_S <- mfbvar_obj$prior_S
  prior_Pi_Omega <- mfbvar_obj$prior_Pi_Omega
  prior_Pi_mean <- mfbvar_obj$prior_Pi_mean
  prior_psi_Omega <- mfbvar_obj$prior_psi_Omega
  prior_psi_mean <- mfbvar_obj$prior_psi_mean

  freq <- mfbvar_obj$mfbvar_prior$freq
  n_q <- sum(freq == "q")
  T_b <- max(which(!apply(apply(Y[, freq == "m"], 2, is.na), 1, any)))
  Lambda <- mfbvar:::build_Lambda(ifelse(freq == "q", "average", freq), n_lags)
  Lambda_ <- mfbvar:::build_Lambda(rep("average", n_q), 3)

  ################################################################
  ### Initialize
  Pi_red <- array(NA, dim = c(n_vars, n_vars * n_lags, n_reps))
  Sigma_red <- array(NA, dim = c(n_vars, n_vars, n_reps))
  Z_red <- array(NA, dim = c(n_T, n_vars, n_reps))

  Pi_red[, , 1] <- post_Pi_mean
  Sigma_red[, , 1] <- post_Sigma
  Z_red[, , 1] <- apply(mfbvar_obj$Z, c(1, 2), mean)

  roots <- vector("numeric", n_reps)
  num_tries <- roots

  ################################################################
  ### Compute terms which do not vary in the sampler

  # Create D (does not vary in the sampler), and find roots of Pi
  D <- mfbvar:::build_DD(d = d, n_lags = n_lags)

  # For the posterior of Pi
  inv_prior_Pi_Omega <- solve(prior_Pi_Omega)
  Omega_Pi <- inv_prior_Pi_Omega %*% prior_Pi_mean

  Z_1 <- Z_red[1:n_lags, , 1]

  mZ <- Y - d %*% t(post_psi)
  mZ <- as.matrix(mZ)
  demeaned_z0 <- Z_1 - d[1:n_lags, ] %*% t(post_psi)
  d_post_psi <- d %*% t(post_psi)
  ################################################################
  ### Reduced Gibbs step
  mod_red <- estimate_mfbvar(
    mfbvar_prior = mfbvar_obj$mfbvar_prior,
    prior = "ss",
    variance = "iw",
    init = list(psi = post_psi),
    fixate = list(psi = TRUE)
  )
  Z_red <- mod_red$Z
  ################################################################
  ### For the likelihood calculation
  mZ <- Y - d %*% t(post_psi)
  mZ <- mZ[-(1:n_lags), ]
  demeaned_z0 <- Z[1:n_lags, , 1] - d[1:n_lags, ] %*% t(post_psi)
  h0 <- matrix(t(demeaned_z0), ncol = 1)
  h0 <- h0[(n_vars * n_lags):1, , drop = FALSE] # have to reverse the order
  Pi_comp <- mfbvar:::build_companion(post_Pi_mean, n_vars = n_vars, n_lags = n_lags)
  Q_comp <- matrix(0, ncol = n_vars * n_lags, nrow = n_vars * n_lags)
  Q_comp[1:n_vars, 1:n_vars] <- t(chol(post_Sigma))
  P0 <- matrix(0, n_lags * n_vars, n_lags * n_vars)

  ################################################################
  ### Final calculations
  lklhd <- sum(c(mfbvar:::loglike(
    Y = as.matrix(mZ), Lambda = Lambda,
    Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_,
    n_vars = n_vars, n_comp = n_lags * n_vars,
    z0 = h0, P0 = P0
  )[-1]))
  eval_prior_Pi_Sigma <- mfbvar:::dnorminvwish(
    X = t(post_Pi_mean), Sigma = post_Sigma,
    M = prior_Pi_mean, P = prior_Pi_Omega,
    S = prior_S, v = n_vars + 2
  )
  eval_prior_psi <- mfbvar:::dmultn(
    x = post_psi, m = prior_psi_mean,
    Sigma = prior_psi_Omega
  )
  eval_log_RB <- mfbvar:::eval_Pi_Sigma_RaoBlack(
    Z_array = Z_red, d = d,
    post_psi_center = t(post_psi),
    post_Pi_center = post_Pi_mean,
    post_Sigma_center = post_Sigma,
    post_nu = post_nu,
    prior_Pi_mean = prior_Pi_mean,
    prior_Pi_Omega = prior_Pi_Omega,
    prior_S = prior_S, n_vars = n_vars,
    n_lags = n_lags, n_reps = n_reps
  )
  const <- median(eval_log_RB)
  eval_RB_Pi_Sigma <- log(mean(exp(eval_log_RB - const))) + const
  eval_marg_psi <- log(mean(mfbvar:::eval_psi_MargPost(
    Pi_array = Pi, Sigma_array = Sigma,
    Z_array = Z,
    post_psi_center = post_psi,
    prior_psi_mean = prior_psi_mean,
    prior_psi_Omega = prior_psi_Omega,
    D_mat = D, n_determ = n_determ,
    n_vars = n_vars, n_lags = n_lags,
    n_reps = n_reps
  )))

  mdd_estimate <- c(lklhd + eval_prior_Pi_Sigma + eval_prior_psi - (eval_RB_Pi_Sigma + eval_marg_psi))

  return(list(lklhd = lklhd, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma, eval_prior_psi = eval_prior_psi, eval_RB_Pi_Sigma = eval_RB_Pi_Sigma, eval_marg_psi = eval_marg_psi, log_mdd = mdd_estimate))
}

#' Estimate marginal data density in Minnesota MF-BVAR
#'
#' This function provides the possibility to estimate the log marginal density (up to a constant) using the Minnesota MF-BVAR.
#' @rdname mdd.minn
#' @param mfbvar_obj An object of class \code{mfbvar} containing the results
#' @param quarterly_cols numeric vector with positions of quarterly variables
#' @param p_trunc \code{1-p_trunc} is the degree of truncation (i.e.,
#' \code{p_trunc = 1} is no truncation) in the truncated normal distribution
#' @keywords internal
#' @noRd
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

  postsim <- vapply(
    1:n_reps, function(x) {
      Z_comp <- mfbvar:::build_Z(z = Z[, , x], n_lags = n_lags)
      XX <- Z_comp[-nrow(Z_comp), ]
      XX <- cbind(XX, 1)
      YY <- Z_comp[-1, 1:n_vars]

      XXt.XX <- crossprod(XX)
      XXt.XX.inv <- chol2inv(chol(XXt.XX))
      Pi_sample <- XXt.XX.inv %*% crossprod(XX, YY)

      # Posterior moments of Pi
      post_Pi_Omega <- chol2inv(chol(inv_prior_Pi_Omega + XXt.XX))
      post_Pi <- post_Pi_Omega %*% (Omega_Pi + crossprod(XX, YY))
      S <- crossprod(YY - XX %*% Pi_sample)
      Pi_diff <- prior_Pi_mean - Pi_sample
      post_S <- prior_S + S + t(Pi_diff) %*% chol2inv(chol(prior_Pi_Omega + XXt.XX.inv)) %*% Pi_diff
      return(mfbvar:::dmatt(YY, XX %*% prior_Pi_mean, chol2inv(chol(diag(nrow(YY)) + XX %*% prior_Pi_Omega %*% t(XX))), prior_S, prior_nu))
    },
    numeric(1)
  )

  temp <- apply(Z[-(1:n_lags), , ], 3, function(x) x[is.na(c(mfbvar_obj$Y[-(1:n_lags), ]))])

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
    return(log_mdd = mean(mdd) + sum(!is.na(Y[-(1:n_lags), mfbvar_obj$freq == "q"])) * log(3))
  }
}
