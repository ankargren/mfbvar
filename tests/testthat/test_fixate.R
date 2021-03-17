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

  test_wrapper <- function(prior_obj, prior, variance, params) {
    set.seed(10)
    mod <- estimate_mfbvar(prior_obj, prior, variance)

    for (i in params) {
      set.seed(10)
      init <- list(mod[[i]][,,10])
      fixate <- list(TRUE)
      names(init) <- i
      names(fixate) <- i
      cat(prior, variance, i, "\n")
      mod2 <- estimate_mfbvar(prior_obj, prior, variance,
                              init = init,
                              fixate = fixate)
      expect_equal(mod[[i]][,,10], mod2[[i]][,,1])
      expect_equal(mod2[[i]][,,1], mod2[[i]][,,10])
    }
  }


  testthat::skip_on_cran()
  # minn iw
  params <- c("Z", "Pi", "Sigma")
  prior <- "minn"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # ss iw
  params <- c("psi", "Pi", "Sigma", "Z")
  prior <- "ss"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # ssng iw
  params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "Sigma")
  prior <- "ssng"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # minn csv
  params <- c("Z", "Pi", "Sigma", "phi", "sigma", "latent", "latent0")
  prior <- "minn"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

  # ss csv
  params <- c("Z", "psi", "Pi", "Sigma", "phi", "sigma", "latent", "latent0")
  prior <- "minn"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

  # ssng csv
  params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "Sigma",
              "phi", "sigma", "latent", "latent0")
  prior <- "ssng"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

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

