library(mfbvar)
context("Plots")
test_that("Forecasts (minn)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_init(
    Y = Y, freq = c(rep("m", 4), "q"),
    n_lags = 4, n_burnin = 10, n_reps = 10
  )
  prior_obj <- set_prior_minn(prior_obj)
  prior_intervals <- matrix(c(
    6, 7,
    0.1, 0.2,
    0, 0.5,
    -0.5, 0.5,
    0.4, 0.6
  ), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
    d = "intercept",
    prior_psi_mean = prior_psi_mean,
    prior_psi_Omega = prior_psi_Omega
  )

  testthat::skip_on_cran()
  set.seed(10)
  mod_minn <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "minn",
    n_fcst = 12
  )

  expect_error(plot(mod_minn), NA)
  expect_error(plot(mod_minn, plot_start = "2013-07-31"), NA)

  rownames(Y) <- as.character(floor_date(as_date(rownames(Y)), unit = "month"))
  set.seed(10)
  mod_minn <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "minn",
    n_fcst = 12,
    Y = Y
  )
  expect_error(plot(mod_minn), NA)
  expect_error(plot(mod_minn, plot_start = "2013-07-01"), NA)

  rownames(Y) <- NULL

  set.seed(10)
  mod_minn <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "minn",
    n_fcst = 12,
    Y = Y
  )
  expect_error(plot(mod_minn))
})

test_that("Forecasts (ss)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_init(
    Y = Y, freq = c(rep("m", 4), "q"),
    n_lags = 4, n_burnin = 10, n_reps = 10
  )
  prior_obj <- set_prior_minn(prior_obj)

  prior_intervals <- matrix(c(
    6, 7,
    0.1, 0.2,
    0, 0.5,
    -0.5, 0.5,
    0.4, 0.6
  ), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
    d = "intercept",
    prior_psi_mean = prior_psi_mean,
    prior_psi_Omega = prior_psi_Omega
  )

  testthat::skip_on_cran()
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss",
    n_fcst = 12
  )

  expect_error(plot(mod_ss), NA)
  expect_error(plot(mod_ss, plot_start = "2013-07-31"), NA)

  rownames(Y) <- as.character(floor_date(as_date(rownames(Y)), unit = "month"))
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss", Y = Y,
    n_fcst = 12
  )
  expect_error(plot(mod_ss), NA)
  expect_error(plot(mod_ss, plot_start = "2013-07-01"), NA)

  rownames(Y) <- NULL

  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss", Y = Y,
    n_fcst = 12
  )
  expect_error(plot(mod_ss))
})

test_that("Prior", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_init(
    Y = Y, freq = c(rep("m", 4), "q"),
    n_lags = 4, n_burnin = 10, n_reps = 10
  )
  expect_error(plot(prior_obj), NA)

  rownames(Y) <- NULL
  prior_obj <- set_init(
    Y = Y, freq = c(rep("m", 4), "q"),
    n_lags = 4, n_burnin = 10, n_reps = 10
  )
  plot(prior_obj)
})

test_that("varplot", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_init(
    Y = Y, freq = c(rep("m", 4), "q"),
    n_lags = 4, n_burnin = 10, n_reps = 10
  )
  prior_obj <- set_prior_minn(prior_obj)
  prior_obj <- set_prior_fsv(prior_obj, n_fac = 1)
  prior_obj <- set_prior_csv(prior_obj)

  prior_intervals <- matrix(c(
    6, 7,
    0.1, 0.2,
    0, 0.5,
    -0.5, 0.5,
    0.4, 0.6
  ), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
    d = "intercept",
    prior_psi_mean = prior_psi_mean,
    prior_psi_Omega = prior_psi_Omega
  )

  testthat::skip_on_cran()
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss",
    variance = "fsv"
  )
  expect_error(varplot(mod_ss, variables = "gdp"), NA)

  rownames(Y) <- NULL
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss", Y = Y,
    variance = "fsv"
  )
  expect_error(varplot(mod_ss, variables = "gdp"), NA)

  colnames(Y) <- NULL
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss", Y = Y,
    variance = "fsv"
  )
  expect_error(varplot(mod_ss, variables = 1), NA)

  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss", Y = Y,
    variance = "csv"
  )
  expect_error(varplot(mod_ss, variables = 1), NA)
})

test_that("Weekly-Monthly plots", {
  set.seed(10237)
  Y <- matrix(rnorm(400), 100, 4)
  Y[setdiff(1:100, seq(4, 100, by = 4)), 4] <- NA

  prior_obj <- set_init(
    Y = Y, freq = c(rep("w", 3), "m"),
    n_lags = 4, n_reps = 10, n_burnin = 10
  )
  prior_obj <- set_prior_minn(prior_obj)

  prior_intervals <- matrix(c(
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5
  ), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
    d = "intercept",
    prior_psi_mean = prior_psi_mean,
    prior_psi_Omega = prior_psi_Omega
  )
  prior_obj <- set_prior_csv(prior_obj)

  testthat::skip_on_cran()
  set.seed(10)
  mod_ss <- estimate_mfbvar(
    mfbvar_prior = prior_obj, prior = "ss",
    variance = "csv"
  )
  expect_error(varplot(mod_ss, variables = 1), NA)
  expect_error(plot(mod_ss))
})
