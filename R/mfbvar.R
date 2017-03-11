#' MF-BVAR with a steady-state prior
#'
#' The main function for running the MF-SS-BVAR.
#'
#' @templateVar Y TRUE
#' @templateVar d TRUE
#' @templateVar d_fcst TRUE
#' @templateVar Lambda TRUE
#' @templateVar prior_Pi_AR1 TRUE
#' @templateVar lambda1 TRUE
#' @templateVar lambda2 TRUE
#' @templateVar prior_nu TRUE
#' @templateVar prior_psi_mean TRUE
#' @templateVar prior_psi_Omega TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_fcst TRUE
#' @templateVar n_burnin TRUE
#' @templateVar n_reps TRUE
#' @template man_template
#' @details \code{mfbvar} calls \code{\link{gibbs_sampler}} (implemented in C++)

mfbvar <- function(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2, prior_nu = NULL, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps) {

  stopifnot(is.matrix(d), is.matrix(prior_psi_Omega))
  stopifnot(is.data.frame(Y) || is.matrix(Y), all(apply(Y, 2, is.numeric)))
  stopifnot(is.vector(n_lags), is.vector(n_burnin), is.vector(n_reps), is.vector(lambda1), is.vector(lambda2), is.vector(n_fcst))
  stopifnot(nrow(Y) == nrow(d), ncol(Y) == length(prior_Pi_AR1), ncol(Y) == length(prior_psi_mean),
            ncol(Y) == sqrt(prod(dim(prior_psi_Omega))))

  if (!is.numeric(Lambda)) {
    if (is.character(Lambda) && length(Lambda) == ncol(Y)) {
      Lambda <- build_Lambda(Lambda, n_lags)
    } else {
      stop("Lambda must be either a matrix of correct size or a character vector with the aggregation use.")
    }
  }

  if (!is.null(n_fcst)) {
    if (nrow(d_fcst) != n_fcst) {
      stop(paste0("d_fcst has ", nrow(d_fcst), " rows, but n_fcst is ", n_fcst, "."))
    } else {
      if (is.null(rownames(d_fcst))) {
        names_fcst <- paste0("fcst_", 1:n_fcst)
      } else {
        names_fcst <- rownames(d_fcst)
      }

    }
  }

  fun_call <- match.call()
  if (is.null(rownames(Y))) {
    names_row <- 1:nrow(Y)
  } else {
    names_row <- rownames(Y)
  }

  if (is.null(colnames(Y))) {
    names_col <- 1:col(Y)
  } else {
    names_col <- colnames(Y)
  }

  if (is.null(colnames(d))) {
    names_determ <- paste0("d", 1:ncol(d))
  } else {
    names_determ <- colnames(d)
  }

  original_Y <- Y
  Y <- as.matrix(Y)
  n_vars <- ncol(Y)
  n_T <- nrow(Y)

  # Set priors
  if (is.null(prior_nu)) {
    prior_nu <- n_vars + 2
  }
  priors <- prior_Pi_Sigma(lambda1 = lambda1, lambda2 = lambda2, prior_Pi_AR1 = prior_Pi_AR1, Y = Y,
                           n_lags = n_lags, prior_nu = prior_nu)
  prior_Pi_mean <- priors$prior_Pi_mean
  prior_Pi_Omega <- priors$prior_Pi_Omega
  prior_S <- priors$prior_S

  # For the smoothing
  cat(paste0("############################################\n   Running the burn-in sampler with ", n_burnin, " draws\n\n"))
  start_burnin <- Sys.time()
  burn_in <-  gibbs_sampler(Y = Y, d = d, d_fcst = NULL, Lambda = Lambda, prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = prior_Pi_Omega,
                            prior_S = prior_S, prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega,
                            n_fcst = NULL, n_reps = n_burnin, init_Pi = t(prior_Pi_mean), init_Sigma = NULL, init_psi = prior_psi_mean,
                            init_Z = NULL, smooth_state = FALSE, check_roots = TRUE)
  end_burnin <- Sys.time()
  time_diff <- end_burnin - start_burnin
  cat(paste0("\n   Time elapsed for drawing ", n_burnin, " times for burn-in: ", signif(time_diff, digits = 1), " ",
             attr(time_diff, "units"), "\n"))
  cat(paste0("\n   Moving on to ",
             n_reps, " replications in the main chain\n", ifelse(!is.null(n_fcst), paste0("   Making forecasts ", n_fcst, " steps ahead"), NULL), "\n\n"))

  main_run <- gibbs_sampler(Y = Y, d = d, d_fcst = d_fcst, Lambda = Lambda, prior_Pi_mean = prior_Pi_mean, prior_Pi_Omega = prior_Pi_Omega,
                prior_S = prior_S, prior_nu = prior_nu, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega,
                n_fcst = n_fcst, n_reps = n_burnin, init_Pi  = burn_in$Pi[,,dim(burn_in$Pi)[3]], init_Sigma = burn_in$Sigma[,,dim(burn_in$Sigma)[3]],
                init_psi = burn_in$psi[dim(burn_in$psi)[1],], init_Z   = burn_in$Z[,,dim(burn_in$Z)[3]], smooth_state = FALSE, check_roots = TRUE)
  main_run$call <- fun_call
  time_diff <- Sys.time() - start_burnin
  cat(paste0("\n  Total time elapsed: ", signif(time_diff, digits = 1), " ",
             attr(time_diff, "units")))

  if (!is.null(n_fcst)) {
    rownames(main_run$Z_fcst)[1:main_run$n_lags] <- names_row[(main_run$n_T-main_run$n_lags+1):main_run$n_T]
    rownames(main_run$Z_fcst)[(main_run$n_lags+1):(main_run$n_fcst+main_run$n_lags)] <- names_fcst
    colnames(main_run$Z_fcst) <- names_col

  }

  main_run$names_row <- names_row
  main_run$names_col <- names_col
  main_run$names_fcst <- names_fcst
  main_run$names_determ <- names_determ


  class(main_run) <- "mfbvar"
  return(main_run)

}

#' Printing method for class mfbvar
#'
#' The method for printing mfbvar objects.
#'
#' @param x object of class mfbvar
#' @param ... Currently not in use.
#' @template man_template

print.mfbvar <- function(x, ...){
  cat(paste0("Mixed-frequency steady-state BVAR with:\n", ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "), "\n", nrow(x$prior_Pi_mean)/ncol(x$prior_Pi_mean), " lags\n",
             nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecasted\n",
             x$n_reps, " draws used in main chain"))
}

#' Plotting method for class mfbvar
#'
#' The method for plotting mfbvar objects.
#'
#' @param x object of class mfbvar
#' @param plot_start Time period (number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param ss_level A vector with the lower and upper quantiles for the posterior steady-state intervals.
#' @param pred_level A vector with the lower and upper quantiles for the forecast intervals.
#' @param ... Currently not in use.
#' @template man_template
plot.mfbvar <- function(x, plot_start = NULL, ss_level = c(0.025, 0.975),
                        pred_level = c(0.10, 0.90), ...){
  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (!is.null(x$n_fcst)) {
      plot_range <- max(x$n_T-x$n_fcst*5, 0):x$n_T
    }  else {
      plot_range <- 1:x$n_T
    }
  } else {
    plot_range <- plot_start:x$n_T
  }

  if (is.null(ss_level)) {
    ss_level <- c(0.025, 0.975)
  }
  if (is.null(pred_level)) {
    pred_level <- c(0.10, 0.90)
  }

  ss_lower  <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
  ss_median <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
  ss_upper  <- x$d[plot_range, ] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))
  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  names_row <- if (is.null(x$names_row)) 1:x$n_T else x$names_row
  ss <- data.frame(expand.grid(time = plot_range, names_col = names_col), lower = c(ss_lower), median = c(ss_median),
                   upper = c(ss_upper))
  ss$value <- c(as.matrix(x$Y[plot_range,]))
  ss <- na.omit(ss)

  p <- ggplot(ss, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "grey90"), alpha =1) +
    geom_line(aes(y = value))

  if (!is.null(x$n_fcst)) {
    preds <- predict(x, pred_quantiles = c(pred_level[1], 0.5, pred_level[2]))
    fcst <- data.frame(expand.grid(time = (x$n_T+1):(x$n_T+x$n_fcst), names_col = names_col),
                       lower = c(preds[[1]]), median = c(preds[[2]]),
                     upper = c(preds[[3]]))
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(fcst, data.frame(time = (1:nrow(x$Y))[last_pos[i]],
                                              names_col = names_col[i],
                                              lower = x$Y[last_pos[i], i],
                                              median = x$Y[last_pos[i], i],
                                              upper  = x$Y[last_pos[i], i]))
    }
    p <- p + geom_ribbon(data = fcst, aes(ymin = lower, ymax = upper,
                                          fill = "#bdbdbd"), alpha = 1) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = "Note: The forecasts are for the underlying variable.")
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }



  p <- p + facet_wrap(~names_col, scales = "free_y") +
    scale_fill_manual(values = c("#bdbdbd", "grey90"),
                      label = c(paste0("Prediction (", 100*(pred_level[2]-pred_level[1]), "%)"),
                                paste0("Steady state (", 100*(ss_level[2]-ss_level[1]), "%)"))) +
    labs(fill = "Intervals",
         title = "Forecasts and steady state intervals",
         y = "Value",
         x = "Time") +
    theme_minimal() +
    theme(legend.position="bottom")
  breaks <- ggplot_build(p)$layout$panel_ranges[[1]]$x.major_source
  breaks <- breaks[-which(!(breaks %in% 1:x$n_T))]
  p + scale_x_continuous(breaks = breaks,
                         labels = names_row[breaks])
}

#' Summary method for class mfbvar
#'
#' The method for summarizing mfbvar objects (currently the same as \code{print.mfbvar}).
#'
#' @param object object of class mfbvar
#' @param ... Currently not in use.
#' @template man_template
#'
summary.mfbvar <- function(object, ...) {
  post_Pi <- apply(object$Pi, c(1, 2), mean)
  rownames(post_Pi) <- object$names_col
  colnames(post_Pi) <- paste0(rep(object$names_col, object$n_lags), ".", rep(1:object$n_lags, each = object$n_vars))
  post_Sigma <- apply(object$Sigma, c(1, 2), mean)
  rownames(post_Sigma) <- object$names_col
  colnames(post_Sigma) <- object$names_col
  post_psi <- matrix(colMeans(object$psi), ncol = object$n_determ)
  rownames(post_psi) <- object$names_col
  colnames(post_psi) <- object$names_determ
  print(object, ...)
  cat("\n\n#########################\nPosterior means computed\n\nPi:\n")
  print(post_Pi)
  cat("\n\n Sigma:\n")
  print(post_Sigma)
  cat("\n\n Psi:\n")
  print(post_psi)
  ret_list <- list(post_Pi = post_Pi, post_Sigma = post_Sigma, post_Psi = post_psi)
}

#' Predict method for class mfbvar
#'
#' The method for predicting mfbvar objects.
#'
#' @param object object of class mfbvar
#' @param pred_quantiles The quantiles of the posterior predictive distribution to use.
#' @param tidy If results should be tidy or not.
#' @param ... Currently not in use.
#' @template man_template
#' @details Note that this requires that forecasts were made in the original \code{mfbvar} call.
#'
predict.mfbvar <- function(object, pred_quantiles = c(0.10, 0.50, 0.90), tidy = FALSE, ...) {
  if (is.null(object$n_fcst)) {
    stop("No forecasts exist in the provided object.")
  }

  if (tidy == FALSE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[-(1:object$n_lags),,-1], c(1, 2), quantile, prob = xx))
    names(ret_list) <- paste0("quantile_", pred_quantiles*100)
    return(ret_list)
  } else if (tidy == TRUE) {
    ret_list <- lapply(pred_quantiles, function(xx) apply(object$Z_fcst[-(1:object$n_lags),,-1], c(1, 2), quantile, prob = xx))
    ret_tidy <- cbind(value = unlist(ret_list), expand.grid(fcst_date = rownames(object$Z_fcst[-(1:object$n_lags),,2]),
                                        variable = object$names_col,
                                        quantile = pred_quantiles))
    return(ret_tidy)
  }

}
