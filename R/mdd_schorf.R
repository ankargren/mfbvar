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
mdd_schorf <- function(monthly_cols, postsim, Z, n_T) {

  n_TT <- (n_T %% 3) + 1
  lstate0 <- Z[n_TT:n_T, monthly_cols, , drop = FALSE]

  lstateA <- c()
  for (i in seq_along(monthly_cols)) {
    lstate1 <- t(lstate0[seq(from = 1, to = n_T-n_TT+1, by = 3), i, ])
    lstate2 <- t(lstate0[seq(from = 2, to = n_T-n_TT+1, by = 3), i, ])
    lstate3 <- t(lstate0[seq(from = 3, to = n_T-n_TT+1, by = 3), i, ])

    # I don't know why we do this!
    temp <- cbind(lstate3 - lstate2, lstate2 - lstate1)
    lstateA <- cbind(lstateA, temp)
  }

  n_para <- ncol(lstateA)
  n_simul <- nrow(lstateA)

  drawmean <- colMeans(lstateA)
  drawsig <- crossprod(lstateA)/n_simul - tcrossprod(drawmean)


  svd_res <- svd(drawsig, nu = ncol(drawsig))
  up <- svd_res$u
  sp <- rbind(diag(svd_res$d), matrix(0, nrow = dim(up)[1]- length(svd_res$d), ncol = length(svd_res$d)))
  vp <- svd_res$v

  sp_inv <- diag(ifelse(diag(sp) > 1e-12, 1/diag(sp), 0))
  drawsiglndet <- sum(log(ifelse(diag(sp) > 1e-12, diag(sp), 1)))

  # for (j in 1:nrow(drawsig)) {
  #   if (sp[j, j] > 1e-12) {
  #     drawsiginv[j, j] <- 1/sp[j, j]
  #     drawsigdet[j, j] <- sp[j, j]
  #   } else {
  #     drawsigdet[j, j] <- 1
  #   }
  # }


  drawsiginv <- vp %*% tcrossprod(sp_inv, up)

  p <- seq(from = 0.1, to = 0.9, by = 0.1)
  pcrit <- qchisq(p, df = n_para)

  paradev  <- lstateA - kronecker(matrix(1, n_simul, 1), matrix(drawmean, nrow = 1))
  quadpara <- rowSums((paradev %*% drawsiginv) * paradev)

  invlike <- matrix(NA, n_simul, length(p))
  for (i in seq_along(p)) {
    for (j in 1:n_simul) {
      lnfpara <- -0.5*n_para*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara[j] - log(p[i])
      indpara <- quadpara[j] < pcrit[i]
      invlike[j, i] <- exp(lnfpara - (postsim[j]-postsim[1])) * indpara
    }
    meaninvlike <- colMeans(invlike)
    mdd <- -log(meaninvlike)
  }
  return(log_mdd = mdd)
}
