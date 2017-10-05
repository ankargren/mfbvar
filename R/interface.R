#' Create a new object of the S3 class MFBVAR (the class name can be changed)
#'
#' An object of the S3 class MFBVAR will be created inside the function.
#'
#' The members (attributes) of the object give the settings for the mfbvar model.
#' The functions in the package take the object as the input.
#' The members of the object can be initialized by the function, and can also be changed by the following "set..." functions.
#'
#'
set_prior <- function(Y, d = NULL, d_fcst = NULL, Lambda, prior_Pi_AR1 = rep(0, ncol(Y)), lambda1 = 0.2, lambda2 = 1, lambda3 = NULL, prior_nu = NULL, prior_psi_mean = NULL, prior_psi_Omega = NULL, n_lags, n_fcst = 0, n_burnin, n_reps, verbose = FALSE) {
  if (!is.matrix(Y)) {
    if (!is.data.frame(Y)) {
      stop(paste0("Y is of class ", class(Y), ", but must be matrix or data frame."))
    } else {
      Y <- as.matrix(Y)
    }
  }

  if (hasArg(d)) {
    if (!is.matrix(d)) {
      stop("d is of class", class(d), ", but must be a matrix.")
    }

    if (nrow(Y) != nrow(d)) {
      stop("Y has ", nrow(Y), "rows and d ", nrow(d), "rows, but they must be equal.")
    }

    if (hasArg(d_fcst)) {
      if (!is.matrix(d_fcst)) {
        stop("d is of class", class(d), ", but must be a matrix.")
      }

      if (ncol(d) != ncol(d_fcst)) {
        stop("d has", ncol(d), " columns and d_fcst", ncol(d_fcst), ".")
      }
    }
  }

  if (!is.matrix(Lambda) && !is.vector(Lambda)) {
    stop("Lambda is of class ", class(Lambda), ", but it must be a matrix or a character vector.")
  }
  if (is.matrix(Lambda)) {
    if (nrow(Lambda) != ncol(Y)) {
      stop("Lambda has ", nrow(Lambda), " and Y has ", ncol(Y), " columns. The dimensions must match.")
    }
    if (!(ncol(Lambda %% ncol(Y)))) {
      stop("The column dimension of Lambda is ", ncol(Lambda), ", but it must be a multiple of the number of variables (", ncol(Y), ").")
    }
  }
  if (is.vector(Lambda) !all(Lambda %in% c("identity", "average", "triangular"))) {
    stop("Valid aggregations are 'identity', 'average' and 'triangular', but specification includes ", Lambda[!which(Lambda %in% c("identity", "average", "triangular"))], ".")
  }



  if (hasArg(prior_Pi_AR1)) {
    if (is.vector(prior_Pi_AR1)) {
      if (length(prior_Pi_AR1) == 1) {
        cat("prior_Pi_AR1: Recycling ", prior_Pi_AR1, " to use as prior mean for AR(1) coefficients for all variables.\n")
        prior_Pi_AR1 <- rep(prior_Pi_AR1, ncol(Y))
      }
    } else {
      stop("prior_Pi_AR1 must be a vector, but is now of class ", class(prior_Pi_AR1), ".")
    }
  } else {
    cat("prior_Pi_AR1: Using 0 as prior mean for AR(1) coefficients for all variables.\n")
    prior_Pi_AR1 <- rep(0, ncol(Y))
  }

  if (hasArg(lambda1)) {
    if (!is.vector(lambda1) || length(lambda1)) {
      stop("lambda1 must be a vector with a single element.")
    }
  } else {
    cat("lambda2: Using the default", lambda1, "as the value for the lag decay hyperparameter.\n")
  }

  if (hasArg(lambda2)) {
    if (!is.vector(lambda2) || length(lambda2)) {
      stop("lambda2 must be a vector with a single element.")
    }
  } else {
    cat("lambda2: Using the default", lambda2, "as the value for the lag decay hyperparameter.\n")
  }

  if (hasArg(lambda3)) {
    if (!is.vector(lambda3) || length(lambda3)) {
      stop("lambda3 must be a vector with a single element.")
    }
  }

  if (hasArg(prior_nu)) {
    if (!is.vector(prior_nu) || length(prior_nu)) {
      stop("prior_nu must be a vector with a single element.")
    }
  }

  if (hasArg(prior_psi_mean)) {
    if (!(is.vector(prior_psi_mean) || is.matrix(prior_psi_mean))) {
      stop("prior_psi_mean must be a vector or matrix with one row or column.")
    }
    if (is.vector(prior_psi_mean)) {
      if (length(prior_psi_mean) != ncol(Y)) {
        stop("prior_psi_mean has ", length(prior_psi_mean), " elements, but there are ", ncol(Y), " variables in Y.")
      }
    }
    if (is.matrix(prior_psi_mean)) {
      if (!any(dim(prior_psi_mean) == 1)) {
        stop("prior_psi_mean must be a vector or matrix with one row or column.")
      } else {
        prior_psi_mean <- c(prior_psi_mean)
      }
    }
  }

  if (hasArg(prior_psi_Omega)) {
    if (!is.matrix(prior_psi_Omega)) {
      stop("prior_psi_Omega must be a matrix.")
    } else {
      if (dim(prior_psi_Omega)[1] != dim(prior_psi_Omega)[2]) {
        stop("prior_psi_Omega must be a positive-definite symmetric matrix.")
      }
    }
  }

  if (hasArg(n_lags)) {
    if (!is.vector(n_lags) || length(n_lags)) {
      stop("n_lags must be a vector with a single element.")
    }
    if (is.vector(Lambda)) {
      Lambda <- build_Lambda(Lambda, n_lags)
    }
  } else {
    stop("n_lags: No lag length specified.\n")
  }

  if (hasArg(n_fcst)) {
    if (!is.vector(n_fcst) || length(n_fcst)) {
      stop("n_fcst must be a vector with a single element.")
    }
  } else {
    cat("n_fcst: Using the default", n_fcst, "for the number of forecasts to compute.\n")
  }

  if (hasArg(n_burnin)) {
    if (!is.vector(n_burnin) || length(n_burnin)) {
      stop("n_burnin must be a vector with a single element.")
    }
  } else {
    stop("n_burnin: Number of burn-in draws to use not specified.\n")
  }

  if (hasArg(n_reps)) {
    if (!is.vector(n_reps) || length(n_reps)) {
      stop("n_reps must be a vector with a single element.")
    }
  } else {
    stop("n_reps: Number of draws to use in main chain not specified.\n")
  }

  ret <- list(Y = Y, d = d, d_fcst = d_fcst, Lambda = Lambda, prior_Pi_AR1 = prior_Pi_AR1, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3,
              prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega, n_lags = n_lags, n_fcst = n_fcst, n_burnin = n_burnin,
              n_reps = n_reps, verbose = verbose)
  class(ret) <- "mfbvar_prior"

  return(ret)
}

#' Set attributes or members in an object of the class MFBVAR
#'
#' @param obj an object of the class MFBVAR
#' @param ... the attributes and their values to be changed
#'
update_prior <- function(obj, ...)
{
  if(class(obj) != "mfbvar_prior") {
    stop("The object to be updated is not of class mfbvar_prior.")
  }

  tmp <- list(...)
  nam <- names(tmp)
  for(iter in 1:length(tmp)) eval(parse(text=paste0("obj$",nam[iter]," = tmp[[iter]]")))

  return(obj)
}


estimate_mfbvar <- function(mfbvar_prior = NULL, prior_type = c("ss", "minn"), ...) {
  if (hasArg(mfbvar_prior)) {
    if (class(mfbvar_prior) != "mfbvar_prior") {
      stop("mfbvar_prior must be of class mfbvar_prior.")
    } else {
      mfbvar_prior <- update_prior(mfbvar_prior, ...)
    }
  } else {
    mfbvar_prior <- set_prior(...)
  }


  if (!(prior_type %in% c("ss", "minn"))) {
    stop("prior_type must be 'ss' or 'minn'.")
  }

  n_vars <- ncol(mfbvar_prior$Y)
  n_T <- nrow(mfbvar_prior$Y)

  if (prior_type == "ss") {
    if (is.null(mfbvar_prior$prior_nu)) {
      cat("prior_nu: Using the default n_vars + 2 prior for prior_nu.\n")
      mfbvar_prior$prior_nu <- n_vars + 2
    }

    if (mfbvar_prior$n_fcst > 0 && nrow(mfbvar_prior$d_fcst) != mfbvar_prior$n_fcst) {
      stop("d_fcst has ", nrow(mfbvar_prior$d_fcst), " rows, but n_fcst is ", mfbvar_prior$n_fcst, ".")
    }

    priors <- prior_Pi_Sigma(lambda1 = mfbvar_prior$lambda1, lambda2 = mfbvar_prior$lambda2, prior_Pi_AR1 = mfbvar_prior$prior_Pi_AR1, Y = mfbvar_prior$Y,
                             n_lags = mfbvar_prior$n_lags, prior_nu = mfbvar_prior$prior_nu)
    prior_Pi_mean <- priors$prior_Pi_mean
    prior_Pi_Omega <- priors$prior_Pi_Omega
    prior_S <- priors$prior_S

    if (mfbvar_prior$n_fcst > 0) {
      if (is.null(rownames(mfbvar_prior$d_fcst))) {
        names_fcst <- paste0("fcst_", 1:mfbvar_prior$n_fcst)
      } else {
        names_fcst <- rownames(mfbvar_prior$d_fcst)
      }
    } else {
      names_fcst <- NULL
    }

    if (is.null(colnames(mfbvar_prior$d))) {
      names_determ <- paste0("d", 1:ncol(mfbvar_prior$d))
    } else {
      names_determ <- colnames(mfbvar_prior$d)
    }

  } else {
    if (is.null(mfbvar_prior$lambda3)) {
      mfbvar_prior$lambda3 <- 10000
      cat("lambda3: Using default value of 10,000.")
    }

    if (mfbvar_prior$n_fcst > 0) {
      names_fcst <- paste0("fcst_", 1:mfbvar_prior$n_fcst)
    } else {
      names_fcst <- NULL
    }
  }

  if (is.null(rownames(mfbvar_prior$Y))) {
    names_row <- 1:nrow(mfbvar_prior$Y)
  } else {
    names_row <- rownames(mfbvar_prior$Y)
  }

  if (is.null(colnames(mfbvar_prior$Y))) {
    names_col <- 1:col(mfbvar_prior$Y)
  } else {
    names_col <- colnames(mfbvar_prior$Y)
  }



  if (mfbvar_prior$verbose) {
    cat(paste0("############################################\n   Running the burn-in sampler with ", mfbvar_prior$n_burnin, " draws\n\n"))
    start_burnin <- Sys.time()
  }
  if (prior_type == "ss") {
    burn_in <-  gibbs_sampler(Y = mfbvar_prior$Y, d = mfbvar_prior$d, d_fcst = NULL, Lambda = mfbvar_prior$Lambda, prior_Pi_mean = prior_Pi_mean,
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S, prior_nu = mfbvar_prior$prior_nu,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega,
                              n_fcst = NULL, n_reps = mfbvar_prior$n_burnin, smooth_state = FALSE, check_roots = TRUE, verbose = mfbvar_prior$verbose)
  } else {
    burn_in <-  gibbs_sampler_schorf(mfbvar_prior$Y, mfbvar_prior$Lambda, mfbvar_prior$prior_Pi_AR1, mfbvar_prior$lambda1, mfbvar_prior$lambda2,
                                     mfbvar_prior$lambda3, mfbvar_prior$n_lags, n_fcst = NULL, mfbvar_prior$n_burnin,
                                     smooth_state = FALSE, check_roots = TRUE, verbose = mfbvar_prior$verbose)
  }

  if (mfbvar_prior$verbose) {
    end_burnin <- Sys.time()
    time_diff <- end_burnin - start_burnin
    cat(paste0("\n   Time elapsed for drawing ", mfbvar_prior$n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
    cat(paste0("\n   Moving on to ",
               mfbvar_prior$n_reps, " replications in the main chain\n", ifelse(mfbvar_prior$n_fcst > 0, paste0("   Making forecasts ", mfbvar_prior$n_fcst, " steps ahead"), ""), "\n\n"))
  }

  if (prior_type == "ss") {
    main_run <- gibbs_sampler(Y = mfbvar_prior$Y, d = mfbvar_prior$d, d_fcst = mfbvar_prior$d_fcst, Lambda = mfbvar_prior$Lambda, prior_Pi_mean = prior_Pi_mean,
                              prior_Pi_Omega = prior_Pi_Omega, prior_S = prior_S, prior_nu = mfbvar_prior$prior_nu,
                              prior_psi_mean = mfbvar_prior$prior_psi_mean, prior_psi_Omega = mfbvar_prior$prior_psi_Omega, n_fcst = mfbvar_prior$n_fcst,
                              n_reps = mfbvar_prior$n_reps, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                              init_psi = burn_in$psi[dim(burn_in$psi)[1],], init_Z = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = FALSE, check_roots = TRUE,
                              verbose = mfbvar_prior$verbose)
  } else {
    main_run <- gibbs_sampler_schorf(mfbvar_prior$Y, mfbvar_prior$Lambda, mfbvar_prior$prior_Pi_AR1, mfbvar_prior$lambda1, mfbvar_prior$lambda2, mfbvar_prior$lambda3,
                                     mfbvar_prior$n_lags, mfbvar_prior$n_fcst, mfbvar_prior$n_reps, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]],
                                     init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]], init_Z = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = FALSE,
                                     check_roots = TRUE, mfbvar_prior$verbose)
  }


  if (mfbvar_prior$verbose) {
    time_diff <- Sys.time() - start_burnin
    cat(paste0("\n   Total time elapsed: ", signif(time_diff, digits = 1), " ",
               attr(time_diff, "units"), "\n"))
  }

  if (mfbvar_prior$n_fcst == 0) {
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(main_run$n_T-main_run$n_lags+1):main_run$n_T]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col
  }


  main_run$names_row <- names_row
  main_run$names_col <- names_col
  main_run$names_fcst <- names_fcst
    main_run$mfbvar_prior <- mfbvar_prior

  dimnames(main_run$Z) <- list(time = names_row,
                               variable = names_col,
                               iteration = 1:mfbvar_prior$n_reps)
  dimnames(main_run$Pi) <- list(dep = names_col,
                                indep = paste0(rep(names_col, mfbvar_prior$n_lags), ".l", rep(1:mfbvar_prior$n_lags, each = n_vars)),
                                iteration = 1:mfbvar_prior$n_reps)
  dimnames(main_run$Sigma) <- list(names_col,
                                   names_col,
                                   iteration = 1:mfbvar_prior$n_reps)

  if (prior_type == "ss") {
    main_run$names_determ <- names_determ
    n_determ <- dim(mfbvar_prior$d)[2]
    dimnames(main_run$psi) <- list(iteration = 1:mfbvar_prior$n_reps,
                                   param = paste0(rep(names_col, n_determ), ".", rep(names_determ, each = n_vars)))
  }

  class(main_run) <- "mfbvar"
  return(main_run)
}

