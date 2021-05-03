library(mfbvar)
context("Predict")
test_that("Forecasts (mf)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  rownames(Y) <- as.character(lubridate::floor_date(lubridate::ymd(rownames(Y)),
                                                    unit = "months"))
  prior_obj <- set_init(Y = Y, freq = c(rep("m", 4), "q"),
                         n_lags = 4, n_burnin = 10, n_reps = 10)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               0.4, 0.6), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj, d = "intercept",
                            prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega)
  prior_obj <- set_prior_minn(prior_obj)

  testthat::skip_on_cran()
  set.seed(10)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj,
                              prior = "minn",
                              n_fcst = 12)
  expect_equal(predict(mod_minn) %>%
    dplyr::filter(variable == "gdp") %>%
    .$median,
  c(median(colMeans(mod_minn$Z_fcst[2:4,5,])),
  median(colMeans(mod_minn$Z_fcst[5:7,5,])),
  median(colMeans(mod_minn$Z_fcst[8:10,5,])),
  median(colMeans(mod_minn$Z_fcst[11:13,5,])),
  median(colMeans(mod_minn$Z_fcst[14:16,5,]))))

  expect_equal(predict(mod_minn, aggregate_fcst = FALSE) %>%
    dplyr::filter(variable == "gdp") %>%
    .$median,
    c(apply(mod_minn$Z_fcst[-1,5,], 1, median), use.names = FALSE))


  Y_list <- c(lapply(Y[,1:4],
                     function(x) ts(x, frequency = 12, start = c(1996, 8))),
              list(gdp = ts(Y[seq(from = 2, to = nrow(Y), by = 3), 5],
                            frequency = 4, start = c(1996, 3))))

  set.seed(10237)
  prior_obj2 <- set_init(Y = Y_list, n_lags = 4, n_burnin = 10, n_reps = 10)
  prior_obj2 <- set_prior_minn(prior_obj2)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               0.4, 0.6), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj2 <- set_prior_ss(prior_obj2, d = "intercept",
                             prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega)
  set.seed(10)
  mod_minn2 <- estimate_mfbvar(mfbvar_prior = prior_obj2, prior = "minn", n_fcst = 12)

  expect_equal(predict(mod_minn), predict(mod_minn2))

})

test_that("Forecasts (monthly)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  rownames(Y) <- as.character(lubridate::floor_date(lubridate::ymd(rownames(Y)),
                                                    unit = "months"))
  prior_obj <- set_init(Y = na.omit(Y[, -5]), freq = rep("m", 4),
                         n_lags = 4, n_burnin = 10, n_reps = 10)
  prior_obj <- set_prior_minn(prior_obj)
  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj, d = "intercept",
                         prior_psi_mean = prior_psi_mean,
                         prior_psi_Omega = prior_psi_Omega)

  testthat::skip_on_cran()
  set.seed(10)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj,
                              prior = "minn",
                              n_fcst = 12)
  fcst <- predict(mod_minn)
  fcst <- fcst[fcst$variable == "eti", ]$median
  expect_equal(fcst,
               c(apply(mod_minn$Z_fcst[-(1:4),4,], 1, median), use.names = FALSE))

  Y_list <- lapply(na.omit(Y[,1:4]),
                   function(x) ts(x, frequency = 12, start = c(1996, 8)))

  set.seed(10237)
  prior_obj2 <- set_init(Y = Y_list, n_lags = 4, n_burnin = 10, n_reps = 10)
  prior_obj2 <- set_prior_minn(prior_obj2)
  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj2 <- set_prior_ss(prior_obj2,
                             d = "intercept",
                             prior_psi_mean = prior_psi_mean,
                             prior_psi_Omega = prior_psi_Omega)
  set.seed(10)
  mod_minn2 <- estimate_mfbvar(mfbvar_prior = prior_obj2,
                               prior = "minn",
                               n_fcst = 12)

  expect_equal(predict(mod_minn), predict(mod_minn2))
})

test_that("Forecasts (quarterly)", {
  set.seed(10237)
  Y <- mfbvar::mf_sweden
  rownames(Y) <- as.character(lubridate::floor_date(lubridate::ymd(rownames(Y)),
                                                    unit = "months"))
  prior_obj <- set_init(Y = na.omit(Y[seq(2, nrow(Y), by = 3), ]),
                        freq = rep("q", 5),
                         n_lags = 4, n_burnin = 10, n_reps = 10)
  prior_obj <- set_prior_minn(prior_obj)

  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               1, 3), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
                            d = "intercept",
                            prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega)

  testthat::skip_on_cran()
  set.seed(10)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj,
                              prior = "minn", n_fcst = 12)
  fcst <- predict(mod_minn)
  fcst <- fcst[fcst$variable == "eti", ]$median
  expect_equal(fcst,
               as.numeric(apply(mod_minn$Z_fcst[-(1:4),4,], 1, median)))

  Y <- na.omit(Y[seq(2, nrow(Y), by = 3), ])
  Y_list <- lapply(Y, function(x) ts(x, frequency = 4, start = c(1996, 3)))

  set.seed(10237)
  prior_obj2 <- set_init(Y = Y_list, n_lags = 4, n_burnin = 10, n_reps = 10)
  prior_obj2 <- set_prior_minn(prior_obj2)
  prior_intervals <- matrix(c( 6,   7,
                               0.1, 0.2,
                               0,   0.5,
                               -0.5, 0.5,
                               1, 3), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj2 <- set_prior_ss(prior_obj2,
                             d = "intercept",
                             prior_psi_mean = prior_psi_mean,
                             prior_psi_Omega = prior_psi_Omega)
  set.seed(10)
  mod_minn2 <- estimate_mfbvar(mfbvar_prior = prior_obj2,
                               prior = "minn",
                               n_fcst = 12)

  expect_equal(predict(mod_minn), predict(mod_minn2))
})

test_that("Forecasts (weekly-monthly)", {
  set.seed(10237)
  Y <- matrix(rnorm(400), 100, 4)
  Y[setdiff(1:100,seq(4, 100, by = 4)), 4] <- NA

  prior_obj <- set_init(Y = Y, freq = c(rep("w", 3), "m"),
                         n_lags = 4, n_reps = 10, n_burnin = 10)
  prior_obj <- set_prior_minn(prior_obj)
  prior_intervals <- matrix(c(
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5,
    -0.5, 0.5), ncol = 2, byrow = TRUE)
  psi_moments <- interval_to_moments(prior_intervals)
  prior_psi_mean <- psi_moments$prior_psi_mean
  prior_psi_Omega <- psi_moments$prior_psi_Omega
  prior_obj <- set_prior_ss(prior_obj,
                            d = "intercept",
                            prior_psi_mean = prior_psi_mean,
                            prior_psi_Omega = prior_psi_Omega)

  testthat::skip_on_cran()
  set.seed(10)
  mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj,
                              prior = "minn",
                              n_fcst = 12)
  fcst <- predict(mod_minn, aggregate_fcst = FALSE)
  fcst <- fcst[fcst$variable == 4, ]$median
  expect_equal(fcst,
               c(apply(mod_minn$Z_fcst[-(1:4),4,], 1, median),
                 use.names = FALSE))

})
