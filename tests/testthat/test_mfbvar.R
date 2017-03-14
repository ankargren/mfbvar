library(mfbvar)
context("Output")
test_that("Output correct", {
  TT <- 200
  n_vars <- 3
  set.seed(100)

  Y <- matrix(0, 2*TT, n_vars)
  Phi <- matrix(c(0.3, 0.1, 0.2, 0.3, 0.3, 0.6, 0.2, 0.2, 0.3), 3, 3)
  for (i in 2:(2*TT)) {
    Y[i, ] <- Phi %*% Y[i-1,] + rnorm(n_vars)
  }
  Y[, n_vars] <- zoo::rollapply(Y[, n_vars], 3, mean, fill = NA, align = "right")
  Y <- Y[-(1:TT),]
  Y[setdiff(1:TT, seq(1, TT, 3)), n_vars] <- NA

  dates <- paste(rep(2000:2017, each = 12), "-", 1:12, sep = "")
  Y <- as.data.frame(Y)
  names_row <- dates[1:nrow(Y)]
  rownames(Y) <- names_row
  names_col <- c("GDP", "Infl", "Interest")
  colnames(Y) <- names_col

  n_burnin <- 200
  n_reps <- 200
  n_fcst <- 8
  n_lags <- 4
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  d <- matrix(1, nrow = n_T, ncol = 1, dimnames = list(1:nrow(Y), "const"))
  d_fcst <- matrix(1, nrow = n_fcst, ncol = 1,
                   dimnames = list(dates[(nrow(Y)+1):(nrow(Y)+n_fcst)], "const"))

  prior_nu <- n_vars + 2
  prior_Pi_AR1 <- c(0, 0, 0)
  lambda1 <- 0.1
  lambda2 <- 1


  prior_psi_mean <- c(0, 0, 0)
  prior_psi_Omega <- c(0.5, 0.5, 0.5)
  prior_psi_Omega <- diag((prior_psi_Omega / (qnorm(0.975, mean = 0, sd = 1)*2))^2)

  Lambda <- build_Lambda(c("identity", "identity", "average"), n_lags)

  set.seed(10237)
  mfbvar_obj <- mfbvar(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2,
                       prior_nu, prior_psi_mean, prior_psi_Omega,
                       n_lags, n_fcst, n_burnin, n_reps, verbose = TRUE)

  expect_equal(unname(c(mfbvar_obj$Y[, 1])), c(mfbvar_obj$Z[, 1, 100]))

  expect_equal(unname(c(mfbvar_obj$Y[, 1])), c(mfbvar_obj$Z[, 1, 100]))

  expect_equal_to_reference(mfbvar_obj$Z[,3, 100], "Z_output.rds")
  expect_equal_to_reference(mfbvar_obj$Pi[,, 100], "Pi_output.rds")
  expect_equal_to_reference(mfbvar_obj$psi[100, ], "psi_output.rds")
  expect_equal_to_reference(mfbvar_obj$Sigma[,, 100], "Sigma_output.rds")


})
