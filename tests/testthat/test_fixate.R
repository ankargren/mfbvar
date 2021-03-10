library(mfbvar)
context("Fixated values")
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
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(Pi = mod$Pi[,,10]),
                          fixate = list(Pi = TRUE))
  expect_equal(mod$Pi[,,10], mod2$Pi[,,1])
  expect_equal(mod2$Pi[,,1], mod2$Pi[,,10])

  # Initial Sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(Sigma = mod$Sigma[,,10]),
                          fixate = list(Sigma = TRUE))
  expect_equal(mod$Sigma[,,10], mod2$Sigma[,,1])
  expect_equal(mod2$Sigma[,,1], mod2$Sigma[,,10])

  # Initial Z
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(Z = mod$Z[,,10]),
                          fixate = list(Z = TRUE))
  expect_equal(mod$Z[,,10], mod2$Z[,,1])
  expect_equal(mod2$Z[,,1], mod2$Z[,,10])

  # SS
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "ss", "iw")

  # Initial steady state
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(psi = mod$psi[,,10]),
                          fixate = list(psi = TRUE))
  expect_equal(mod$psi[,,10], mod2$psi[,,1])
  expect_equal(mod2$psi[,,1], mod2$psi[,,10])

  # SSNG
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "ssng", "iw")

  # Initial omega
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(omega = mod$omega[,,10]),
                          fixate = list(omega = TRUE))
  expect_equal(mod$omega[,,10], mod2$omega[,,1])
  expect_equal(mod2$omega[,,1], mod2$omega[,,10])

  # Initial phi_mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(phi_mu = mod$phi_mu[,,10]),
                          fixate = list(phi_mu = TRUE))
  expect_equal(mod$phi_mu[,,10], mod2$phi_mu[,,1])
  expect_equal(mod2$phi_mu[,,1], mod2$phi_mu[,,10])

  # Initial lambda_mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(lambda_mu = mod$lambda_mu[,,10]),
                          fixate = list(lambda_mu = TRUE))
  expect_equal(mod$lambda_mu[,,10], mod2$lambda_mu[,,1])
  expect_equal(mod2$lambda_mu[,,1], mod2$lambda_mu[,,10])

  # CSV
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "minn", "csv")

  # Initial phi
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(phi = mod$phi[,,10]),
                          fixate = list(phi = TRUE))
  expect_equal(mod$phi[,,10], mod2$phi[,,1])
  expect_equal(mod2$phi[,,1], mod2$phi[,,10])

  # Initial sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(sigma = mod$sigma[,,10]),
                          fixate = list(sigma = TRUE))
  expect_equal(mod$sigma[,,10], mod2$sigma[,,1])
  expect_equal(mod2$sigma[,,1], mod2$sigma[,,10])

  # Initial latent
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(latent = mod$latent[,,10]),
                          fixate = list(latent = TRUE))
  expect_equal(mod$latent[,,10], mod2$latent[,,1])
  expect_equal(mod2$latent[,,1], mod2$latent[,,10])

  # FSV
  set.seed(10)
  mod <- estimate_mfbvar(prior_obj, "minn", "fsv")

  # Initial mu
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(mu = mod$mu[,,10]),
                          fixate = list(mu = TRUE))
  expect_equal(mod$mu[,,10], mod2$mu[,,1])
  expect_equal(mod2$mu[,,1], mod2$mu[,,10])

  # Initial sigma
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(sigma = mod$sigma[,,10]),
                          fixate = list(sigma = TRUE))
  expect_equal(mod$sigma[,,10], mod2$sigma[,,1])
  expect_equal(mod2$sigma[,,1], mod2$sigma[,,10])

  # Initial phi
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(phi = mod$phi[,,10]),
                          fixate = list(phi = TRUE))
  expect_equal(mod$phi[,,10], mod2$phi[,,1])
  expect_equal(mod2$phi[,,1], mod2$phi[,,10])

  # Initial f
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(f = mod$f[,,10]),
                          fixate = list(f = TRUE))
  expect_equal(mod$f[,,10], mod2$f[,,1])
  expect_equal(mod2$f[,,1], mod2$f[,,10])

  # Initial latent0
  set.seed(10)
  mod2 <- estimate_mfbvar(prior_obj, "minn", "iw",
                          init = list(latent0 = mod$latent0[,,10]),
                          fixate = list(latent0 = TRUE))
  expect_equal(mod$latent0[,,10], mod2$latent0[,,1])
  expect_equal(mod2$latent0[,,1], mod2$latent0[,,10])
})

