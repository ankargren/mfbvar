#' Estimate MDD in Schorfheide-Song MFBVAR
#'
#' This function provides the possibility to estimate the log marginal density (up to a constant) using the Schorfheide-Song MFBVAR.
#'
#' @templateVar monthly_cols TRUE
#' @templateVar postsim TRUE
#' @templateVar Z TRUE
#' @templateVar n_T TRUE
#' @template man_template
#' @return The log marginal data density estimate (bar a constant)
mdd_schorf <- function(mfbvar_obj, monthly_cols) {
  postsim <- mfbvar_obj$lnpYY[-1]
  Z <- mfbvar_obj$Z
  n_T <- mfbvar_obj$n_T
  n_reps <- mfbvar_obj$n_reps -1

  temp <- apply(Z[-(1:mfbvar_obj$n_lags), , -1], 3, function(x) x[is.na(c(mfbvar_obj$Y[-(1:mfbvar_obj$n_lags),]))])
  n_para <- nrow(temp)

  drawmean <- matrix(rowMeans(temp), ncol = 1)
  drawsig <- tcrossprod(temp/sqrt(ncol(temp))) - tcrossprod(drawmean)
  drawsiginv <- chol2inv(chol(drawsig))
  drawsiglndet <- log(det(drawsiginv))
  # n_TT <- (n_T %% 3) + 1
  # lstate0 <- Z[n_TT:n_T, monthly_cols, , drop = FALSE]
  #
  # lstateA <- c()
  # for (i in seq_along(monthly_cols)) {
  #   lstate1 <- t(lstate0[seq(from = 1, to = n_T-n_TT+1, by = 3), i, ])
  #   lstate2 <- t(lstate0[seq(from = 2, to = n_T-n_TT+1, by = 3), i, ])
  #   lstate3 <- t(lstate0[seq(from = 3, to = n_T-n_TT+1, by = 3), i, ])
  #
  #   # I don't know why we do this!
  #   temp <- cbind(lstate3 - lstate2, lstate2 - lstate1)
  #   lstateA <- cbind(lstateA, temp)
  # }
  #
  #   n_para <- ncol(lstateA)
  #   n_simul <- nrow(lstateA)
  #
  #   drawmean <- colMeans(lstateA)
  #   drawsig <- crossprod(lstateA)/n_simul - tcrossprod(drawmean)
  #
  #
  #   svd_res <- svd(drawsig, nu = ncol(drawsig))
  #   up <- svd_res$u
  #   sp <- rbind(diag(svd_res$d), matrix(0, nrow = dim(up)[1]- length(svd_res$d), ncol = length(svd_res$d)))
  #   vp <- svd_res$v
  #
  #   sp_inv <- diag(ifelse(diag(sp) > 1e-12, 1/diag(sp), 0))
  #   drawsiglndet <- sum(log(ifelse(diag(sp) > 1e-12, diag(sp), 1)))
  #
  #   # for (j in 1:nrow(drawsig)) {
  #   #   if (sp[j, j] > 1e-12) {
  #   #     drawsiginv[j, j] <- 1/sp[j, j]
  #   #     drawsigdet[j, j] <- sp[j, j]
  #   #   } else {
  #   #     drawsigdet[j, j] <- 1
  #   #   }
  #   # }
  #
  #
  #   drawsiginv <- vp %*% tcrossprod(sp_inv, up)
  #
  #   p <- seq(from = 0.1, to = 0.9, by = 0.1)
  #   pcrit <- qchisq(p, df = n_para)
  #
  #   paradev  <- lstateA - kronecker(matrix(1, n_simul, 1), matrix(drawmean, nrow = 1))
  #   quadpara <- rowSums((paradev %*% drawsiginv) * paradev)


  p <- seq(from = 0.1, to = 0.9, by = 0.1)
  pcrit <- qchisq(p, df = nrow(drawmean))

  paradev  <- temp - kronecker(matrix(1, 1, n_reps), drawmean)
  quadpara <- rowSums((t(paradev) %*% drawsiginv) * t(paradev))

  invlike <- matrix(NA, n_reps, length(p))
  for (i in seq_along(p)) {
    for (j in 1:n_reps) {
      lnfpara <- -0.5*n_para*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara[j] - log(p[i])
      indpara <- quadpara[j] < pcrit[i]
      invlike[j, i] <- exp(lnfpara - (postsim[j]-postsim[1])) * indpara
    }
    meaninvlike <- colMeans(invlike)
    mdd <- -log(meaninvlike)
  }
  return(log_mdd = mdd)
}
