parallel_wrapper <- function(data_list, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, n_cores, seed, cluster_type = "PSOCK") {

  # data_list is a list of lists with top components: data (with data) and fcst_date (date at which the fcst is made)
  # The components of data are data frames/matrices with rownames containing dates
  d_list <- lapply(data_list$mf, function(x) matrix(1, nrow = nrow(x), ncol = 1, dimnames = list(time = rownames(x), const = "const")))
  d_fcst_list <- lapply(as.Date(sapply(data_list$mf, function(x) rownames(x)[nrow(x)]), origin = "1970-01-01"),
                        function(x) matrix(1, nrow = n_fcst, ncol = 1, dimnames = list(time = as.character(x + lubridate::days(1) + months(1:n_fcst) - lubridate::days(1)), const = "const")))


  if (cluster_type == "MPI") {
    cl <- parallel::makeCluster(n_cores, type = "MPI")
  } else {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  }

  parallel::clusterSetRNGStream(cl, iseed = seed)

  # Load the package
  parallel::clusterEvalQ(cl, {
    library(mfbvar)
  })
  res <- parallel::clusterMap(cl, fun = worker_fun, Y = data_list$data, d = d_list, d_fcst = d_fcst_list,
                     MoreArgs = list(Lambda = Lambda, prior_Pi_AR1 = prior_Pi_AR1, lambda1_grid = lambda1_grid, lambda2_grid = lambda2_grid,
                                     prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega, n_burnin = n_burnin, n_reps = n_reps))
  parallel::stopCluster(cl)
  return(res)
#
}

worker_fun <- function(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, prior_nu = NULL, prior_psi_mean, prior_psi_Omega, n_burnin, n_reps) {
  n_lags <- ncol(Lambda)/nrow(Lambda)
  n_fcst <- nrow(d_fcst)
  lambda_mat <- expand.grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid)

  mapply_fun <- function(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2, prior_nu, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps) {
    mfbvar_obj <- tryCatch({
      mfbvar(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2, prior_nu, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE)
    }, error = function(cond) {
      NULL
      }
    )
    if (!is.null(mfbvar_obj$Z_fcst)) {
      fcst <- mfbvar_obj$Z_fcst
      mdd_est <- tryCatch({
        mdd1(mfbvar_obj)$log_mdd
      }, error = function(cond) {
        NA
      })
    } else {
      mdd_est <- NA
      fcst <- NULL
    }
    return(list(lambda1 = lambda1, lambda2 = lambda2, log_mdd = mdd_est, fcst = fcst))
  }
  ret <- mapply(mapply_fun, lambda1 = lambda_mat[, 1], lambda2 = lambda_mat[, 2],
                MoreArgs = list(Y = Y, d = d, d_fcst = d_fcst, Lambda = Lambda, prior_Pi_AR1 = prior_Pi_AR1,
                                prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega,
                                n_lags = n_lags, n_fcst = n_fcst, n_burnin = n_burnin, n_reps = n_reps))
}
