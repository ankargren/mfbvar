library(mfbvar)
context("Initial values")
test_that("IW", {
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
                            prior_psi_Omega = prior_psi_Omega, n_fcst = 4)

  testthat::skip_on_cran()
  # MINN
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "minn", "iw")

  # Initial Pi
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw", init = list(Pi = mod$Pi[,,10]))
  expect_false(isTRUE(all.equal(mod$Pi, mod2$Pi)))

  # Initial Sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw", init = list(Sigma = mod$Sigma[,,10]))
  expect_false(isTRUE(all.equal(mod$Pi, mod2$Pi)))

  # Initial Z -- sampled first so shouldn't matter
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw", init = list(Z = mod$Z[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # SS
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "ss", "iw")

  # Initial steady state
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "ss", "iw", init = list(psi = mod$psi[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # SSNG
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "ssng", "iw")

  # Initial omega
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "ssng", "iw", init = list(omega = mod$omega[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial phi_mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "ssng", "iw", init = list(phi_mu = mod$phi_mu[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial phi_mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "ssng", "iw", init = list(lambda_mu = mod$lambda_mu[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # CSV
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "minn", "csv")

  # Initial phi
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "csv", init = list(phi = mod$phi[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "csv", init = list(sigma = mod$sigma[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "csv", init = list(latent = mod$latent[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # FSV
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "minn", "fsv")

  # Initial mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "fsv", init = list(mu = mod$mu[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "fsv", init = list(sigma = mod$sigma[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial phi
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "fsv", init = list(phi = mod$phi[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial f
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "fsv", init = list(f = mod$f[,,10]))
  expect_equal(mod$Pi, mod2$Pi)

  # Initial latent0
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "fsv", init = list(latent0 = mod$latent0[,,10]))
  expect_equal(mod$Pi, mod2$Pi)
})

