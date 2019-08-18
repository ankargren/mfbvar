library(mfbvar)
context("Output")
test_that("Output correct", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000)

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
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", n_fcst = 4)
  mod_ss <- estimate_mfbvar(prior_obj2, "ss")

  mdd_minn <- mdd(mod_minn, p_trunc = 0.5)
  mdd_ss1 <- mdd(mod_ss)

})
context("Prior checks")
test_that("Prior checks correct", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden

  # If Y is not matrix/df
  expect_error(prior_obj <- set_prior(Y = "test", freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))

  # Including d
  expect_error(prior_obj <- set_prior(Y = Y, d = matrix(1, nrow = nrow(Y)-1, 1), freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000))

  # freq
  expect_error(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "s"),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_error(prior_obj <- set_prior(Y = Y, freq = list(c(rep("m", 4), "s")),
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))


  expect_error(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                      aggregation = "triangular",
                                      n_lags = 4, n_burnin = 100, n_reps = 1000))
  expect_error(prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                      aggregation = "average",
                                      n_lags = 2, n_burnin = 100, n_reps = 1000))
  # Using update
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                                        n_lags = 4, n_burnin = 100, n_reps = 1000)
  prior_obj2 <- update_prior(prior_obj, d = "intercept", Y = Y[1:100, ], n_fcst = 4)
  expect_is(prior_obj2$d_fcst, "matrix")
  expect_is(prior_obj2$d, "matrix")

  prior_obj2 <- update_prior(prior_obj, d = "intercept", Y = Y[1:90, ])
  prior_obj2 <- update_prior(prior_obj2, n_fcst = 4)
  expect_is(prior_obj2$d_fcst, "matrix")
  expect_true(all(dim(prior_obj2$d_fcst) == c(4, 1)))
})
