#' Estimate MDD in Schorfheide-Song MFBVAR
#'
#' This function provides the possibility to estimate the log marginal density (up to a constant) using the Schorfheide-Song MFBVAR.
#'
#' @templateVar mfbvar_obj TRUE
#' @template man_template
#' @return The log marginal data density estimate (bar a constant)
mdd_schorf <- function(mfbvar_obj, quarterly_cols, type = "full", ...) {
  postsim <- mfbvar_obj$lnpYY[-1]
  Z <- mfbvar_obj$Z[,, -1]
  n_T <- dim(Z)[1]
  n_reps <- dim(Z)[3]

  temp <- apply(Z[-(1:mfbvar_obj$n_lags), , ], 3, function(x) x[is.na(c(mfbvar_obj$Y[-(1:mfbvar_obj$n_lags),]))])

  if (length(temp) == 0) {
    return(log_mdd = mean(postsim))
  } else {
    if (type == "diff") {
      n_TT <- (n_T %% 3) + 1
      lstate0 <- Z[n_TT:n_T, quarterly_cols, , drop = FALSE]

      lstateA <- c()
      for (i in seq_along(quarterly_cols)) {
        lstate1 <- t(lstate0[seq(from = 1, to = n_T-n_TT+1, by = 3), i, ])
        lstate2 <- t(lstate0[seq(from = 2, to = n_T-n_TT+1, by = 3), i, ])
        lstate3 <- t(lstate0[seq(from = 3, to = n_T-n_TT+1, by = 3), i, ])

        # I don't know why we do this!
        temp <- cbind(lstate3 - lstate2, lstate2 - lstate1)
        lstateA <- cbind(lstateA, temp)
      }

      n_para <- ncol(lstateA)
      n_simul <- nrow(lstateA)

      drawmean <- matrix(colMeans(lstateA), ncol = 1)
      drawsig <- crossprod(lstateA)/n_simul - tcrossprod(drawmean)


      svd_res <- svd(drawsig, nu = ncol(drawsig))
      up <- svd_res$u
      sp <- rbind(diag(svd_res$d), matrix(0, nrow = dim(up)[1]- length(svd_res$d), ncol = length(svd_res$d)))
      vp <- svd_res$v

      sp_inv <- diag(ifelse(diag(sp) > 1e-12, 1/diag(sp), 0))
      drawsiglndet <- sum(log(ifelse(diag(sp) > 1e-12, diag(sp), 1)))

      drawsiginv <- 0*drawsig
      drawsigdet <- 0*drawsig
      for (j in 1:nrow(drawsig)) {
        if (sp[j, j] > 1e-12) {
          drawsiginv[j, j] <- 1/sp[j, j]
          drawsigdet[j, j] <- sp[j, j]
        } else {
          drawsigdet[j, j] <- 1
        }
      }


      drawsiginv <- vp %*% tcrossprod(sp_inv, up)

      paradev  <- lstateA - kronecker(matrix(1, n_simul, 1), matrix(drawmean, nrow = 1))
      quadpara <- rowSums((paradev %*% drawsiginv) * paradev)
    } else if (type == "full") {
      n_para <- nrow(temp)

      drawmean <- matrix(rowMeans(temp), ncol = 1)
      drawsig <- tcrossprod(temp/sqrt(ncol(temp))) - tcrossprod(drawmean)
      drawsiginv <- chol2inv(chol(drawsig))
      drawsiglndet <- as.numeric(determinant(drawsiginv, logarithm = TRUE)$modulus)

      paradev  <- temp - kronecker(matrix(1, 1, n_reps), drawmean)
      quadpara <- rowSums((t(paradev) %*% drawsiginv) * t(paradev))
    }

    p <- seq(from = 0.1, to = 0.9, by = 0.1)
    pcrit <- qchisq(p, df = nrow(drawmean))

    invlike <- matrix(NA, n_reps, length(p))
    indpara <- invlike
    lnfpara <- indpara
    densfac <- -0.5*n_para*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara[1] - log(p) - postsim[1]
    densfac <- -mean(densfac)
    for (i in seq_along(p)) {
      for (j in 1:n_reps) {
        lnfpara[j, i] <- -0.5*n_para*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara[j] - log(p[i])
        indpara[j, i] <- quadpara[j] < pcrit[i]
        invlike[j, i] <- exp(lnfpara[j, i] - postsim[j]+densfac) * indpara[j, i]
      }
      meaninvlike <- colMeans(invlike)
      mdd <- -log(meaninvlike)
    }
    return(log_mdd = mean(mdd))
  }
}
