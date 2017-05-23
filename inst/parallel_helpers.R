parallel_wrapper <- function(data_list, prior, freq, pars, n_cores, cluster_type = "PSOCK", seed, same_seed = FALSE) {


  # data_list is a list of lists with top components: data (with data) and fcst_date (date at which the fcst is made)
  # The components of data are data frames/matrices with rownames containing dates
  d_list <- lapply(data_list$data, function(x) matrix(1, nrow = nrow(x), ncol = 1, dimnames = list(time = rownames(x), const = "const")))
  d_fcst_list <- lapply(as.Date(sapply(data_list$data, function(x) rownames(x)[nrow(x)]), origin = "1970-01-01"),
                          function(x) matrix(1, nrow = n_fcst, ncol = 1, dimnames = list(time = as.character(x + lubridate::days(1) + months(1:n_fcst, abbreviate = TRUE) - lubridate::days(1)), const = "const")))

  if (n_cores > 1) {
    if (cluster_type == "MPI") {
      cl <- parallel::makeCluster(n_cores, type = "MPI")
    } else {
      cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    }


    parallel::clusterExport(cl, varlist = c("same_seed", "seed"))

    if (same_seed == FALSE) {
      parallel::clusterSetRNGStream(cl, iseed = seed)
    }

    # Load the package
    parallel::clusterEvalQ(cl, {
      library(mfbvar)
    })

    res <- parallel::clusterMap(cl, fun = worker_fun, Y = data_list$data, d = d_list, d_fcst = d_fcst_list,
                                MoreArgs = list(prior = prior, freq = freq, pars = pars, n_burnin = n_burnin, n_reps = n_reps,
                                                seed = seed, same_seed = same_seed))
    parallel::stopCluster(cl)
  } else {
    res <- list()
    for (i in seq_along(data_list$data)) {
      print(i)
      res[[i]] <- worker_fun(Y = data_list$data[[i]], d = d_list[[i]], d_fcst = d_fcst_list[[i]], prior = prior, freq = freq, pars = pars,
                           n_burnin = n_burnin, n_reps = n_reps,
                                    seed = seed, same_seed = same_seed)
    }

  }

  return(res)
#
}

worker_fun <- function(Y, prior, freq, d, d_fcst, pars, n_burnin, n_reps, seed, same_seed) {
  n_fcst <- pars$n_fcst
  n_lags <- pars$n_lags
  lambda_mat <- expand.grid(lambda1 = pars$lambda1_grid, lambda2 = pars$lambda2_grid)

  mapply_fun <- function(Y, prior, freq, d, d_fcst, lambda1, lambda2, pars, n_burnin, n_reps, seed, same_seed) {
    for (i in names(pars)) {
      assign(i, pars[[i]])
    }
    if (same_seed) {
      set.seed(seed)
    }
    if (prior == "SS") {
      if (freq == "MF") {
        mfbvar_obj <- #tryCatch({
          mfbvar(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2, prior_nu, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
        #}, error = function(cond) {NULL})
      }
      if (freq == "QF") {
        mfbvar_obj <- #tryCatch({
          qfbvar(Y, d, d_fcst, prior_Pi_AR1, lambda1, lambda2, prior_nu, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
        #}, error = function(cond) {NULL})
      }
    } else {
      if (freq == "MF") {
        mfbvar_obj <- #tryCatch({
          mfbvar_schorf(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
        #}, error = function(cond) {NULL})
      }
      if (freq == "QF") {
        mfbvar_obj <- #tryCatch({
          qfbvar_schorf(Y, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
        #}, error = function(cond) {NULL})
      }
    }

    if (!is.null(mfbvar_obj$Z_fcst)) {
      fcst <- t(mfbvar_obj$Z_fcst[-(1:3),5,])
      if (prior == "SS") {
        if (freq == "MF") {
          mdd_est <- tryCatch({
            mdd1(mfbvar_obj)$log_mdd
          }, error = function(cond) {
            NA
          })
        }
        if (freq == "QF") {
          mdd_est <- tryCatch({
            mdd1_qf(mfbvar_obj)$log_mdd
          }, error = function(cond) {
            NA
          })
        }

      }
      if (prior == "Minn") {
        if (freq == "MF") {
          mdd_est <- tryCatch({
            mdd_schorf(mfbvar_obj)
          }, error = function(cond) {
            NA
          })
        }
        if (freq == "QF") {
          mdd_est <- mfbvar_obj$lnpYY
        }
      }
    }
    else {
      mdd_est <- NA
      fcst <- NULL
    }
    return(list(lambda1 = lambda1, lambda2 = lambda2, log_mdd = c(mdd_est), fcst = fcst))
  }
  ret <- mapply(mapply_fun, lambda1 = lambda_mat[, 1], lambda2 = lambda_mat[, 2],
                MoreArgs = list(Y = Y, prior = prior, freq = freq, d = d, d_fcst = d_fcst, pars = pars, n_burnin = n_burnin, n_reps = n_reps, seed = seed, same_seed = same_seed))
  return(ret)
}


parallel_wrapper_schorf <- function(data_list, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, lambda3, n_lags, n_fcst, n_burnin, n_reps, n_cores, cluster_type = "PSOCK", seed, same_seed = FALSE) {

  if (cluster_type == "MPI") {
    cl <- parallel::makeCluster(n_cores, type = "MPI")
  } else {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  }

  parallel::clusterExport(cl, varlist = c("same_seed", "seed"))

  if (same_seed == FALSE) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  # Load the package
  parallel::clusterEvalQ(cl, {
    library(mfbvar)
  })
  res <- parallel::clusterMap(cl, fun = worker_fun_schorf, Y = data_list$data, MoreArgs = list(Lambda = Lambda, prior_Pi_AR1 = prior_Pi_AR1, lambda1_grid = lambda1_grid,
                                                                                        lambda2_grid = lambda2_grid, lambda3 = lambda3, n_fcst = n_fcst, n_burnin = n_burnin,
                                                                                        n_reps = n_reps, seed = seed, same_seed = same_seed))
  parallel::stopCluster(cl)
  return(res)
  #
}

worker_fun_schorf <- function(Y, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, lambda3, n_fcst, n_burnin, n_reps, seed, same_seed) {
  n_lags <- ncol(Lambda)/nrow(Lambda)
  lambda_mat <- expand.grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid)

  mapply_fun <- function(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_burnin, n_reps, seed, same_seed) {
    if (same_seed) {
      set.seed(seed)
    }
    mfbvar_obj <- tryCatch({
      mfbvar_schorf(Y, Lambda, prior_Pi_AR1, lambda1, lambda2, lambda3, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
    }, error = function(cond) {
      NULL
    }
    )
    if (!is.null(mfbvar_obj$Z_fcst)) {
      fcst <- t(mfbvar_obj$Z_fcst[-(1:3),5,])
      mdd_est <- tryCatch({
        mdd_schorf(mfbvar_obj)
      }, error = function(cond) {
        NA
      })
    } else {
      mdd_est <- NA
      fcst <- NULL
    }
    return(list(lambda1 = lambda1, lambda2 = lambda2, log_mdd = c(mdd_est), fcst = fcst))
  }
  ret <- mapply(mapply_fun, lambda1 = lambda_mat[, 1], lambda2 = lambda_mat[, 2],
                MoreArgs = list(Y = Y, Lambda = Lambda, prior_Pi_AR1 = prior_Pi_AR1, lambda3 = lambda3,
                                n_lags = n_lags, n_fcst = n_fcst, n_burnin = n_burnin, n_reps = n_reps, seed = seed, same_seed = same_seed))
  return(ret)
}

