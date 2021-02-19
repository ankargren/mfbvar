library(mfbvar)
context("MCMC Running")
test_that("Mixed", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                         n_lags = 4, n_burnin = 10, n_reps = 10)

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
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")
})

test_that("Quarterly", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y[seq(2, nrow(Y)-1, by = 3), ], freq = rep("q", 5),
                         n_lags = 4, n_burnin = 10, n_reps = 10)

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
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl", variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")
})

test_that("Monthly", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y[, -5], freq = rep("m", 4),
                         n_lags = 4, n_burnin = 10, n_reps = 10)

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
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl", variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")
})

test_that("Block exogenous 1", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                         n_lags = 4, n_burnin = 10, n_reps = 10,
                         block_exo = 2:3)

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
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv", a = 1/10)
})


test_that("Block exogenous 2", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
                         n_lags = 4, n_burnin = 10, n_reps = 10,
                         block_exo = c("infl", "ip"))

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
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv", a = 1/10)
})

test_that("Weekly-Monthly MCMC", {
  set.seed(10237)
  Y <- matrix(rnorm(400), 100, 4)
  Y[setdiff(1:100,seq(4, 100, by = 4)), 4] <- NA

  prior_obj <- set_prior(Y = Y, freq = c(rep("w", 3), "m"),
                         n_lags = 4, n_reps = 10)

  prior_intervals <- matrix(c(
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4,
                            n_fac = 1)

  testthat::skip_on_cran()
  set.seed(10)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "iw"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "iw"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "iw"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "diffuse"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "diffuse"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "diffuse")

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "csv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "csv"), NA)
  #expect_error(mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "csv"))

  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ss",   variance = "fsv"), NA)
  expect_error(estimate_mfbvar(mfbvar_prior = prior_obj, prior = "ssng", variance = "fsv"), NA)
  #mod <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "dl",   variance = "fsv")
})
