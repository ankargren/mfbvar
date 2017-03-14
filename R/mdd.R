
#' Marginal data density
#'
#' Estimate the marginal data density.
#' @details \code{mdd1} uses method 1, \code{mdd2} uses method 2.
#' @templateVar mfbvar_obj TRUE
#' @templateVar p_trunc TRUE
#' @template man_template
#' @return
#' \code{mdd} returns a list with components being \code{n_reps}-long vectors. These can be used to estimate the MDD.
#' \item{eval_posterior_Pi_Sigma}{Posterior of Pi and Sigma.}
#' \item{data_likelihood}{The likelihood.}
#' \item{eval_prior_Pi_Sigma}{Prior of Pi and Sigma.}
#' \item{eval_prior_psi}{Prior of psi.}
#' \item{psi_truncated}{The truncated psi pdf.}
#' \item{log_mdd}{The mdd estimate (in log).}

mdd2 <- function(mfbvar_obj, p_trunc) {
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
  prior_nu <- mfbvar_obj$prior_nu
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
  eval_prior_Pi_Sigma <- vector("numeric", n_reps)
  eval_prior_psi <- vector("numeric", n_reps)
  psi_truncated <- vector("numeric", n_reps)
  for (r in 1:n_reps) {
    # Demean z, create Z (companion form version)
    demeaned_z <- Z[,, r] - d %*% t(matrix(psi[r, ], nrow = n_vars))
    demeaned_Z <- build_Z(z = demeaned_z, n_lags = n_lags)
    XX <- demeaned_Z[-nrow(demeaned_Z), ]
    YY <- demeaned_Z[-1, 1:n_vars]
    Pi_sample <- solve(crossprod(XX)) %*% crossprod(XX, YY)
    ################################################################
    ### Pi and Sigma step

    # Posterior moments of Pi
    post_Pi_Omega_i <- solve(solve(prior_Pi_Omega) + crossprod(XX))
    post_Pi_i       <- post_Pi_Omega_i %*% (solve(prior_Pi_Omega) %*% prior_Pi_mean + crossprod(XX, YY))

    # Then Sigma
    s_sample  <- crossprod(YY - XX %*% Pi_sample)
    Pi_diff <- prior_Pi_mean - Pi_sample
    post_s_i <- prior_S + s_sample + t(Pi_diff) %*% solve(post_Pi_Omega_i + solve(crossprod(XX))) %*% Pi_diff

    # Set the variables which vary in the Kalman filtering
    mZ <- Y - d %*% t(matrix(psi[r, ], nrow = n_vars))
    mZ <- mZ[-(1:n_lags), ]
    demeaned_z0 <- Z[1:n_lags,, 1] - d[1:n_lags, ] %*% t(matrix(psi[r, ], nrow = n_vars))
    h0 <- matrix(t(demeaned_z0), ncol = 1)
    h0 <- h0[(n_vars*n_lags):1,,drop = FALSE] # have to reverse the order


    eval_posterior_Pi_Sigma[r] <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = post_Pi_i, P = post_Pi_Omega_i, S = post_s_i, v = post_nu)
    data_likelihood[r] <- exp(sum(c(loglike(Y = as.matrix(mZ), Lambda = Lambda, Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_, n_vars = n_vars, n_comp = n_lags * n_vars, z0 = h0, P0 = P0)[-1])))
    eval_prior_Pi_Sigma[r] <- dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = prior_Pi_mean, P = prior_Pi_Omega, S = prior_S, v = prior_nu)
    eval_prior_psi[r] <- dmultn(x = psi[r, ], m = prior_psi_mean, Sigma = prior_psi_Omega)
    psi_truncated[r] <- dnorm_trunc(psi[r, ], post_psi, solve(post_psi_Omega), n_determ*n_vars, p_trunc, chisq_val)

  }

  mdd_num <- eval_posterior_Pi_Sigma * psi_truncated
  log_mdd_den <- log(data_likelihood) + log(eval_prior_Pi_Sigma) + log(eval_prior_psi) # pi_sigma_prior constant? likelihood 0?
  log_mdd <- -log(mean(mdd_num/exp(log_mdd_den)))

  return(list(eval_posterior_Pi_Sigma = eval_posterior_Pi_Sigma, data_likelihood = data_likelihood, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma,
              eval_prior_psi = eval_prior_psi, psi_truncated = psi_truncated, log_mdd = log_mdd))
}

#' @rdname mdd2
#' @return
#' \code{mdd1} returns a list with components (all are currently in logarithms):
#' \item{lklhd}{The likelihood.}
#' \item{eval_prior_Pi_Sigma}{The evaluated prior.}
#' \item{eval_prior_psi}{The evaluated prior of psi.}
#' \item{eval_RB_Pi_Sigma}{The Rao-Blackwellized estimate of the conditional posterior of Pi and Sigma.}
#' \item{eval_marg_psi}{The evaluated marginal posterior of psi.}
#' \item{log_mdd}{The mdd estimate (in log).}
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

  ################################################################
  ### Reduced Gibbs step
  for (r in 2:n_reps) {
    ################################################################
    ### Pi and Sigma step
    #                             (Z_r1,                 d,     psi_r1,            prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots, n_vars, n_lags, n_T)
    Pi_Sigma <- posterior_Pi_Sigma(Z_r1 = Z_red[,, r-1], d = d, psi_r1 = post_psi, prior_Pi_mean, prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S, prior_nu, check_roots = TRUE, n_vars, n_lags, n_T)
    Pi_red[,,r]      <- Pi_Sigma$Pi_r
    Sigma_red[,,r]   <- Pi_Sigma$Sigma_r
    num_tries[r] <- Pi_Sigma$num_try
    roots[r]     <- Pi_Sigma$root

    ################################################################
    ### Smoothing step
    #(Y, d, Pi_r,            Sigma_r,               psi_r,                          Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state)
    Z_res <- posterior_Z(Y, d, Pi_r = Pi_red[,, r], Sigma_r = Sigma_red[,, r], psi_r = post_psi, Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state = FALSE)
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
  lklhd          <- sum(c(loglike(Y = as.matrix(mZ), Lambda = Lambda, Pi_comp = Pi_comp, Q_comp = Q_comp, n_T = n_T_, n_vars = n_vars, n_comp = n_lags * n_vars, z0 = h0, P0 = P0)[-1]))
  eval_prior_Pi_Sigma <- log(dnorminvwish(X = t(post_Pi_mean), Sigma = post_Sigma, M = prior_Pi_mean, P = prior_Pi_Omega, S = prior_S, v = prior_nu))
  eval_prior_psi      <- log(dmultn(x = post_psi, m = prior_psi_mean, Sigma = prior_psi_Omega))
  eval_RB_Pi_Sigma    <- log(mean(eval_Pi_Sigma_RaoBlack(Z_red, d, post_psi, post_Pi_mean, post_Sigma, post_nu, prior_Pi_mean, prior_Pi_Omega, prior_S, n_vars, n_lags, n_reps)))
  eval_marg_psi   <- log(mean(eval_psi_MargPost(Pi, Sigma, Z, post_psi, prior_psi_mean, prior_psi_Omega, D, n_determ, n_vars, n_lags, n_reps)))

  mdd_estimate <- lklhd + eval_prior_Pi_Sigma + eval_prior_psi - (eval_RB_Pi_Sigma + eval_marg_psi)

  return(list(lklhd = lklhd, eval_prior_Pi_Sigma = eval_prior_Pi_Sigma, eval_prior_psi = eval_prior_psi, eval_RB_Pi_Sigma = eval_RB_Pi_Sigma, eval_marg_psi = eval_marg_psi, log_mdd = mdd_estimate))
}

#' Marginal data density grid search
#'
#' Use a grid search to estimate the marginal data density for multiple values, possibly using parallel computation.
#' @details Only one of \code{mfbvar_obj} and the entire suite of \code{Y}, \code{d}, \code{Lambda}, \code{prior_Pi_AR1}, \code{prior_nu}, \code{prior_psi_mean}, \code{prior_psi_Omega}, \code{n_lags}, \code{n_burnin} and \code{n_reps} must be provided. If \code{mfbvar_obj} is not \code{NULL}, its components are extracted.\cr \cr In order to use \code{n_cores > 1}, the \code{parallel} package must be installed.
#' @templateVar mfbvar_obj TRUE
#' @templateVar lambda1_grid TRUE
#' @templateVar lambda2_grid TRUE
#' @templateVar method TRUE
#' @templateVar n_cores TRUE
#' @templateVar p_trunc TRUE
#' @param ... Additional parameters. This includes \code{p_trunc} if \code{method == 2} and if \code{mfbvar_obj} is omitted, it must include all the necessary variables to the \code{mfbvar()} call.
#' @template man_template
#' @return mdd_res A matrix with the results.
#'
mdd_grid <- function(mfbvar_obj = NULL, lambda1_grid, lambda2_grid, method, n_cores = 1, p_trunc = NULL, ...) {
  if (!requireNamespace(c("parallel"), quietly = TRUE) && n_cores > 1) {
    n_cores <- 1
    warnings("The parallel and doParallel are not available. Setting n_cores to 1. To use parallel processing, please install parallel and doParallel.")
  }
  stopifnot((method == 1) || (method == 2))
  if (method == 2) {
    if (is.null(p_trunc)) {
      stop("When using method 2, p_trunc must be supplied.")
    }
  }
  if (is.null(mfbvar_obj)) {
    stopifnot(hasArg(Y), hasArg(d), hasArg(Lambda), hasArg(prior_Pi_AR1),
              hasArg(prior_nu), hasArg(prior_psi_mean), hasArg(prior_psi_Omega),
              hasArg(n_lags), hasArg(n_burnin), hasArg(n_reps))
  } else {
    stopifnot(class(mfbvar_obj) == "mfbvar")
    Y <- mfbvar_obj$Y
    d <- mfbvar_obj$d
    Lambda <- mfbvar_obj$Lambda
    prior_Pi_AR1 <- mfbvar_obj$prior_Pi_AR1
    prior_nu <- mfbvar_obj$prior_nu
    prior_psi_mean <- mfbvar_obj$prior_psi_mean
    prior_psi_Omega <- mfbvar_obj$prior_psi_Omega
    n_lags <- mfbvar_obj$n_lags
    n_burnin <- mfbvar_obj$n_burnin
    n_reps <- mfbvar_obj$n_reps
  }

  mdd_search <- function(lambda1, lambda2, method, p_trunc = NULL) {
    mfbvar_obj <- mfbvar(Y, d, d_fcst = NULL, Lambda, prior_Pi_AR1, lambda1, lambda2,
                         prior_nu, prior_psi_mean, prior_psi_Omega,
                         n_lags, n_fcst = NULL, n_burnin, n_reps, verbose = FALSE)
    if (method == 1) {
      log_mdd <- mdd1(mfbvar_obj)$log_mdd
    } else {
      log_mdd <- mdd2(mfbvar_obj, p_trunc)$log_mdd
    }
    return(log_mdd)
  }

  par_func <- function(j, lambda_mat, method, p_trunc = NULL) {
    lambda1 <- lambda_mat[j, 1]
    lambda2 <- lambda_mat[j, 2]
    cat("Combination", j, "using lambda1 =", lambda1, "lambda2 =", lambda2)
    mdd_res <- tryCatch({mdd_search(lambda1, lambda2, method, p_trunc)},
                        error = function(cond) return(NA),
                        warning = function(cond) return(NA))
    return(c(mdd_res, lambda1, lambda2))
  }

  lambda_mat <- expand.grid(lambda1_grid, lambda2_grid)

  #library(parallel)
  #library(doParallel)
  if (n_cores > 1) {
    cat("Computing log marginal data density\n\nInitiating parallel processing using", n_cores, "cores\n\n")
    cl <- parallel::makePSOCKcluster(n_cores)
    #lusterSetRNGStream(cl, iseed = 2983468)
    parallel::clusterExport(cl = cl, varlist = setdiff(ls(), c("cl", "lambda1", "lambda2")), envir = environment())
    # Load the package
    parallel::clusterEvalQ(cl, {
      library(mfbvar)
    })

    if (method == 1) {
      mdd_res <- parallel::parSapply(cl = cl, X = 1:nrow(lambda_mat), FUN = par_func, lambda_mat = lambda_mat, method = method)
    }
    if (method == 2) {
      mdd_res <- parallel::parSapply(cl = cl, X = 1:nrow(lambda_mat), FUN = par_func, lambda_mat = lambda_mat, method = method, p_trunc = p_trunc)
    }

    parallel::stopCluster(cl)
  } else {
    cat("Computing log marginal data density\n\n")
    if (method == 1) {
      mdd_res <- sapply(X = 1:nrow(lambda_mat), FUN = par_func, lambda_mat = lambda_mat, method = method)
    }
    if (method == 2) {
      mdd_res <- sapply(X = 1:nrow(lambda_mat), FUN = par_func, lambda_mat = lambda_mat, method = method, p_trunc = p_trunc)
    }
  }

  rownames(mdd_res) <- c("log_mdd", "lambda1", "lambda2")
  mdd_ret <- list(mdd_res = mdd_res)
  class(mdd_ret) <- "mdd"
  return(mdd_ret)
}

#' Printing method for class mdd
#'
#' Method for printing mdd objects.
#'
#' @param x object of class mdd
#' @param ... Currently not in use.
print.mdd <- function(x, ...) {
  print(x[["mdd_res"]])
}

#' Plotting method for class mdd
#'
#' Method for plotting mdd objects.
#'
#' @param x object of class mdd
#' @param ... Currently not in use.
plot.mdd <- function(x, ...) {
  plot_df <- data.frame(t(x[["mdd_res"]]))
  lambda1 <- lambda2 <- log_mdd <- NULL
  ggplot(plot_df, aes(x = lambda1, y = lambda2, fill = log_mdd)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradient(low = "grey20", high = "grey90", name = "log(mdd)") +
    coord_fixed(ratio = 1/(c(range(unique(plot_df$lambda2)) %*% c(-1, 1))/c(range(unique(plot_df$lambda1)) %*% c(-1, 1)))) +
    labs(title = "Tile plot of log marginal data density", x = expression(lambda[1]), y = expression(lambda[2]))
}

#' Summary method for class mdd
#'
#' Method for summarizing mdd objects.
#'
#' @param object object of class mdd
#' @param ... Currently not in use.
summary.mdd <- function(object, ...) {
  pos <- which.max(object[["mdd_res"]][1,])
  cat("Highest log marginal data density:", object[["mdd_res"]][1, pos],
      "\nObtained using lambda1 =", object[["mdd_res"]][2, pos], "and lambda2 =", object[["mdd_res"]][3, pos])
}


