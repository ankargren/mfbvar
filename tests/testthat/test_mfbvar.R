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

  mdd_minn <- mdd(mod_minn)
  mdd_ss1 <- mdd(mod_ss, method = 1)
  mdd_ss2 <- mdd(mod_ss, method = 2, p_trunc = 0.5)

  expect_equal(c(mod_ss$Y[!is.na(Y[,1]), 1]), c(mod_ss$Z[!is.na(Y[,1]), 1, 100]))
  expect_equal(c(mod_minn$Y[!is.na(Y[,1]), 1]), c(mod_minn$Z[!is.na(Y[,1]), 1, 100]))

  expect_equal_to_reference(mod_minn$Z[,3, 100], "Z_minn.rds")
  expect_equal_to_reference(mod_minn$Pi[,, 100], "Pi_minn.rds")
  expect_equal_to_reference(mod_minn$Sigma[,, 100], "Sigma_minn.rds")

  expect_equal_to_reference(mod_ss$Z[,3, 100], "Z_ss.rds")
  expect_equal_to_reference(mod_ss$Pi[,, 100], "Pi_ss.rds")
  expect_equal_to_reference(mod_ss$psi[100, ], "psi_ss.rds")
  expect_equal_to_reference(mod_ss$Sigma[,, 100], "Sigma_ss.rds")

  expect_equal_to_reference(mdd_minn, "mdd_minn.rds")
  expect_equal_to_reference(mdd_ss1, "mdd_ss1.rds")
  expect_equal_to_reference(mdd_ss2, "mdd_ss2.rds")
})
