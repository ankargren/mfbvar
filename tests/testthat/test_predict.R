library(mfbvar)
context("Predict")
test_that("Forecasts (mf)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"),
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
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4)

  testthat::skip_on_cran()
  set.seed(100)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", n_fcst = 12)
  expect_equal(predict(mod_minn) %>%
    filter(variable == "gdp") %>%
    pull(median),
  c(median(colMeans(mod_minn$Z_fcst[2:4,5,])),
  median(colMeans(mod_minn$Z_fcst[5:7,5,])),
  median(colMeans(mod_minn$Z_fcst[8:10,5,])),
  median(colMeans(mod_minn$Z_fcst[11:13,5,])),
  median(colMeans(mod_minn$Z_fcst[14:16,5,]))))
})

test_that("Forecasts (monthly)", {
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
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4)

  testthat::skip_on_cran()
  set.seed(100)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", n_fcst = 12)
  expect_equal(predict(mod_minn) %>%
                 filter(variable == "eti") %>%
                 pull(median),
               as.numeric(apply(mod_minn$Z_fcst[-(1:4),4,], 1, median)))
})

test_that("Forecasts (quarterly)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  prior_obj <- set_prior(Y = Y[seq(2, nrow(Y), by = 3), ], freq = rep("q", 5),
                         n_lags = 4, n_burnin = 100, n_reps = 100)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               1, 3), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4)

  testthat::skip_on_cran()
  set.seed(100)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior = "minn", n_fcst = 12)
  expect_equal(predict(mod_minn) %>%
                 filter(variable == "eti") %>%
                 pull(median),
               as.numeric(apply(mod_minn$Z_fcst[-(1:4),4,], 1, median)))
})
