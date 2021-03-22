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
  params <- mfbvar:::get_params_info(minn = TRUE, iw = TRUE)$params
  prior <- "minn"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # ss iw
  params <- mfbvar:::get_params_info(ss = TRUE, iw = TRUE)$params
  prior <- "ss"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # ssng iw
  params <- mfbvar:::get_params_info(ssng = TRUE, iw = TRUE)$params
  prior <- "ssng"
  variance <- "iw"
  test_wrapper(prior_obj, prior, variance, params)

  # minn csv
  params <- mfbvar:::get_params_info(minn = TRUE, csv = TRUE)$params
  prior <- "minn"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

  # ss csv
  params <- mfbvar:::get_params_info(ss = TRUE, csv = TRUE)$params
  prior <- "ss"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

  # ssng csv
  params <- mfbvar:::get_params_info(ss = TRUE, ssng = TRUE, csv = TRUE)$params
  prior <- "ssng"
  variance <- "csv"
  test_wrapper(prior_obj, prior, variance, params)

  # minn diffuse
  params <- mfbvar:::get_params_info(minn = TRUE, diffuse = TRUE)$params
  prior <- "minn"
  variance <- "diffuse"
  test_wrapper(prior_obj, prior, variance, params)

  # ss diffuse
  params <- mfbvar:::get_params_info(ss = TRUE, diffuse = TRUE)$params
  prior <- "ss"
  variance <- "diffuse"
  test_wrapper(prior_obj, prior, variance, params)

  # ssng diffuse
  params <- mfbvar:::get_params_info(ss = TRUE, ssng = TRUE, diffuse = TRUE)$params
  prior <- "ssng"
  variance <- "diffuse"
  test_wrapper(prior_obj, prior, variance, params)

  prior_obj <- update_prior(prior_obj, n_fac = 1)

  test_wrapper_fsv <- function(prior_obj, prior, params) {
    set.seed(10)
    variance <- "fsv"
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
      # Only set-wise fixation is implemented
      if (i %in% c("f", "facload", "latent", "mu", "phi", "sigma")) {
        expect_false(isTRUE(all.equal(mod[[i]][,,10], mod2[[i]][,,1])))
        expect_false(isTRUE(all.equal(mod2[[i]][,,1], mod2[[i]][,,10])))
      } else {
        expect_equal(mod[[i]][,,10], mod2[[i]][,,1])
        expect_equal(mod2[[i]][,,1], mod2[[i]][,,10])
      }
    }
  }

  # minn fsv
  params <- mfbvar:::get_params_info(minn = TRUE, fsv = TRUE)$params
  prior <- "minn"
  test_wrapper_fsv(prior_obj, prior, params)

  # ss fsv
  params <- mfbvar:::get_params_info(ss = TRUE, fsv = TRUE)$params
  prior <- "ss"
  test_wrapper_fsv(prior_obj, prior, params)

  # ssng fsv
  params <- mfbvar:::get_params_info(ss = TRUE, ssng = TRUE, fsv = TRUE)$params
  prior <- "ssng"
  test_wrapper_fsv(prior_obj, prior, params)


})

