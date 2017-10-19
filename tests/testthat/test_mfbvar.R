library(mfbvar)
context("Output")
test_that("Output correct", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  expect_warning(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               0.4, 0.6), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj2 <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                             prior_psi_Omega = prior_psi_Omega, n_fcst = 4)

  expect_true(!is.null(prior_obj2$d_fcst))

  testthat::skip_on_cran()
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior_type = "minn", n_fcst = 4)
  mod_ss <- estimate_mfbvar(prior_obj2, "ss")

  mdd_minn <- mdd(mod_minn, p_trunc = 0.5)
  mdd_ss1 <- mdd(mod_ss, method = 1)
  mdd_ss2 <- mdd(mod_ss, method = 2, p_trunc = 0.5)

  expect_equal(c(mod_ss$Y[!is.na(Y[,1]), 1]), c(mod_ss$Z[!is.na(Y[,1]), 1, 100]))
  expect_equal(c(mod_minn$Y[!is.na(Y[,1]), 1]), c(mod_minn$Z[!is.na(Y[,1]), 1, 100]))

  # saveRDS(mod_minn$Z[,5, 100], file = "tests/testthat/Z_minn.rds")
  # saveRDS(mod_minn$Pi[,, 100], file = "tests/testthat/Pi_minn.rds")
  # saveRDS(mod_minn$Sigma[,, 100], file = "tests/testthat/Sigma_minn.rds")
  #
  # saveRDS(mod_ss$Z[,5, 100], file = "tests/testthat/Z_ss.rds")
  # saveRDS(mod_ss$Pi[,, 100], file = "tests/testthat/Pi_ss.rds")
  # saveRDS(mod_ss$Sigma[,, 100], file = "tests/testthat/Sigma_ss.rds")
  # saveRDS(mod_ss$psi[100,], file = "tests/testthat/psi_ss.rds")
  #
  # saveRDS(mdd_minn, file = "tests/testthat/mdd_minn.rds")
  # saveRDS(mdd_ss1, file = "tests/testthat/mdd_ss1.rds")
  # saveRDS(mdd_ss2, file = "tests/testthat/mdd_ss2.rds")

  expect_equal_to_reference(mod_minn$Z[,5, 100], "Z_minn.rds")
  expect_equal_to_reference(mod_minn$Pi[,, 100], "Pi_minn.rds")
  expect_equal_to_reference(mod_minn$Sigma[,, 100], "Sigma_minn.rds")

  expect_equal_to_reference(mod_ss$Z[,5, 100], "Z_ss.rds")
  expect_equal_to_reference(mod_ss$Pi[,, 100], "Pi_ss.rds")
  expect_equal_to_reference(mod_ss$psi[100, ], "psi_ss.rds")
  expect_equal_to_reference(mod_ss$Sigma[,, 100], "Sigma_ss.rds")

  expect_equal_to_reference(mdd_minn, "mdd_minn.rds")
  expect_equal_to_reference(mdd_ss1, "mdd_ss1.rds")
  expect_equal_to_reference(mdd_ss2, "mdd_ss2.rds")
})
context("Prior checks")
test_that("Prior checks correct", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden

  # If Y is not matrix/df
  expect_error(prior_obj <- set_prior(Y = "test", freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))
  # Still a matrix
  expect_warning(prior_obj <- set_prior(Y = as.ts(Y), freq = c(rep("m", 4), "q"),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  # Including d
  expect_warning(prior_obj <- set_prior(Y = Y, d = "intercept", freq = c(rep("m", 4), "q"),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_warning(prior_obj <- set_prior(Y = Y, d = matrix(1, nrow = nrow(Y), 1), freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_warning(prior_obj <- set_prior(Y = Y, d = cbind(1, 1:nrow(Y)), freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))

  # freq
  expect_error(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "s"),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_error(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4)),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_error(prior_obj <- set_prior(Y = Y, freq = list(c(rep("m", 4), "s")),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))

  # Using update
  expect_warning(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))
  prior_obj <- update_prior(prior_obj, d = "intercept", Y = Y[1:100, ], n_fcst = 4)
  expect_is(prior_obj$d_fcst, "matrix")
  expect_is(prior_obj$d, "matrix")
})
