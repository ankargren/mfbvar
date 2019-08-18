library(mfbvar)
context("MCMC Running")
test_that("Mixed", {
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
  prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4, n_fac = 1)

  testthat::skip_on_cran()
  set.seed(100)
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "iw"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")
})

test_that("Quarterly", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y[seq(2, nrow(Y), by = 3), ], freq = rep("q", 5),
                         n_lags = 4, n_burnin = 100, n_reps = 100)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               0.4, 0.6), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4, n_fac = 1)

  testthat::skip_on_cran()
  set.seed(100)
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl", variance = "iw"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")
})

test_that("Monthly", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y[, -5], freq = rep("m", 4),
                         n_lags = 4, n_burnin = 100, n_reps = 100)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4, n_fac = 1)

  testthat::skip_on_cran()
  set.seed(100)
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl", variance = "iw"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv")
  expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")

  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse")
  mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")
})
