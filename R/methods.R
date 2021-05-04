

#' Print method for mfbvar_prior
#'
#' Printing method for object of class mfbvar_prior, checking if information
#' in the prior is sufficient for estimating models.
#' @param x prior object (class \code{mfbvar_prior})
#' @param ... additional arguments (currently unused)
#' @details The print method checks whether the steady-state and Minnesota
#'   priors can be used with the current specification. This check is minimal in
#'   the sense that it checks only prior elements with no defaults, and it only
#'   checks for estimation and not forecasting (for which the steady-state prior
#'   requires additional information).
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{estimate_mfbvar}}, \code{\link{summary.mfbvar_prior}}
#' @return  No return value, called for side effects.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' print(prior_obj)
print.mfbvar_prior <- function(x, ...) {
  cat("The following elements of the prior object have not been set: \n", names(sapply(x, is.null))[sapply(x, is.null)])
  cat("\n\n")
  cat("Checking if the steady-state prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$d) && !is.null(x$prior_psi_mean) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "d", "prior_psi_mean", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n")
  }

  cat("Checking if a Minnesota-style prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps)) {
    cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps")]
    cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
  }

  # cat("Checking if the Dirichlet-Laplace prior can be used... ")
  if (!is.null(x$Y) && !is.null(x$n_lags) && !is.null(x$n_burnin) && !is.null(x$n_reps) && !is.null(x$a)) {
    # cat("TRUE\n\n")
  } else {
    test_all <- sapply(x, is.null)
    test_sub <- test_all[c("Y", "n_lags", "n_burnin", "n_reps", "a")]
    # cat("FALSE\n Missing elements:", names(test_sub)[which(test_sub)], "\n\n")
  }

  cat("Checking if common stochastic volatility can be used... ")
  if (length(x$prior_phi) == 2 && length(x$prior_sigma2) == 2) {
    cat("TRUE\n\n")
  } else {
    switch(paste0(as.numeric(is.null(x$prior_phi)), as.numeric(is.null(x$prior_sigma2))),
      "01" = cat("FALSE\n Missing element: prior_sigma2 \n\n"),
      "10" = cat("FALSE\n Missing element: prior_phi \n\n"),
      "00" = cat("FALSE\n Missing elements: prior_phi, prior_sigma2 \n\n")
    )
  }

  cat("Checking if factor stochastic volatility can be used... ")
  if (!is.null(x$n_fac)) {
    cat("TRUE\n\n")
  } else {
    cat("FALSE\n Missing element: n_fac \n\n")
  }

  cat("\n")
}

#' Summary method for mfbvar_prior
#'
#' summary method for object of class mfbvar_prior, showing some basic
#' information regarding the contents of the prior.
#' @param object prior object (class \code{mfbvar_prior})
#' @param ... additional arguments (currently unused)
#' @seealso \code{\link{set_prior}}, \code{\link{update_prior}}, \code{\link{estimate_mfbvar}}, \code{\link{print.mfbvar_prior}}
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 100)
#' summary(prior_obj)
summary.mfbvar_prior <- function(object, ...) {
  cat("PRIOR SUMMARY\n")
  cat("----------------------------\n")
  cat("General specification:\n")
  cat("  Y:", ncol(object$Y), "variables,", nrow(object$Y), "time points\n")
  cat("  aggregation:", object$aggregation, "\n")
  freq_count <- vapply(object$freqs, function(x, freq) sum(x == freq), numeric(1), freq = object$freq)
  freqs <- object$freqs
  freqs <- replace(freqs, freqs == "w", "weekly")
  freqs <- replace(freqs, freqs == "m", "monthly")
  freqs <- replace(freqs, freqs == "q", "quarterly")
  if (length(freq_count) == 1) {
    freq_cat <- sprintf("  freq: %d %s variables\n", freq_count, freqs)
  } else {
    freq_cat <- sprintf("  freq: %s variables\n", paste(sprintf("%d %s", rev(freq_count), rev(freqs)), collapse = ", "))
  }
  cat(freq_cat)
  if (length(object$prior_Pi_AR1) <= 5) {
    disp_prior_Pi_AR1 <- object$prior_Pi_AR1
  } else {
    if (length(unique(object$prior_Pi_AR1)) == 1) {
      disp_prior_Pi_AR1 <- object$prior_Pi_AR1[1]
    } else {
      disp_prior_Pi_AR1 <- sprintf("vector with %d elements", length(object$prior_Pi_AR1))
    }
  }
  cat("  prior_Pi_AR1:", disp_prior_Pi_AR1, "\n")
  cat("  lambda1:", object$lambda1, "\n")
  cat("  lambda2:", object$lambda2, "\n")
  cat("  lambda3:", object$lambda3, "\n")
  cat("  lambda4:", ifelse(is.null(object$lambda4), "<missing> (will rely on default)", object$lambda4), "\n")
  cat("  block_exo:", ifelse(is.null(object$block_exo), "0", length(object$block_exo)), "block exogenous variables\n")
  cat("  n_lags:", object$n_lags, "\n")
  cat("  n_fcst:", object$n_fcst, "\n")
  cat("  n_burnin:", object$n_burnin, "\n")
  cat("  n_reps:", object$n_reps, "\n")
  cat("----------------------------\n")
  cat("Steady-state prior:\n")
  cat("  d:", ifelse(is.null(object$d), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(ncol(object$d), "deterministic variables"))), "\n")
  cat("  d_fcst:", ifelse(is.null(object$d_fcst), "<missing>", ifelse(object$intercept_flag, "intercept", paste0(nrow(object$d_fcst), "forecasts, ", ncol(object$d), "deterministic variables"))), "\n")
  cat("  prior_psi_mean:", ifelse(is.null(object$prior_psi_mean), "<missing>", "prior mean vector of steady states"), "\n")
  cat("  prior_psi_Omega:", ifelse(is.null(object$prior_psi_Omega), "<missing>", "prior covariance matrix of steady states"), "\n")
  cat("  check_roots:", object$check_roots, "\n")
  cat("----------------------------\n")
  cat("Hierarchical steady-state prior:\n")
  cat("  s:", ifelse(is.null(object$s), "<missing>", object$s), "\n")
  cat("  c0:", ifelse(is.null(object$prior_ng), "<missing>", object$prior_ng[1]), "\n")
  cat("  c1:", ifelse(is.null(object$prior_ng), "<missing>", object$prior_ng[2]), "\n")
  cat("----------------------------\n")
  # cat("Dirichlet-Laplace prior:\n")
  # cat("  a:", ifelse(is.null(object[["a"]]), "<missing>", object[["a"]]), "\n")
  # cat("----------------------------\n")
  cat("Common stochastic volatility:\n")
  cat(sprintf("  prior_phi: mean = %g, var = %g", object$prior_phi[1], object$prior_phi[2]), "\n")
  cat(sprintf("  prior_sigma2: mean = %g, df = %d", object$prior_sigma2[1], object$prior_sigma2[2]), "\n")
  cat("----------------------------\n")
  cat("Factor stochastic volatility:\n")
  cat("  n_fac:", ifelse(is.null(object$n_fac), "<missing>", object$n_fac), "\n")
  cat("  n_cores:", ifelse(is.null(object$n_cores), "<missing>", object$n_cores), "\n")
  if ("priormu" %in% object$supplied_args) {
    cat("  priormu:", object$priormu, "\n")
  }
  if ("priorphiidi" %in% object$supplied_args) {
    cat("  priorphiidi:", object$priorphiidi, "\n")
  }
  if ("priorphifac" %in% object$supplied_args) {
    cat("  priorphifac:", object$priorphifac, "\n")
  }
  if ("priorsigmaidi" %in% object$supplied_args) {
    if (length(object$priorsigmaidi) == 1) {
      cat("  priorsigmaidi:", object$priorsigmaidi, "\n")
    } else {
      cat("  priorsigmaidi: vector with", length(object$priorsigmaidi), "elements \n")
    }
  }
  if ("priorsigmafac" %in% object$supplied_args) {
    if (length(object$priorsigmafac) == 1) {
      cat("  priorsigmafac:", object$priorsigmafac, "\n")
    } else {
      cat("  priorsigmafac: vector with", length(object$priorsigmafac), "elements \n")
    }
  }
  if ("priorfacload" %in% object$supplied_args) {
    if (length(object$priorfacload) == 1) {
      cat("  priorfacload:", object$priorfacload, "\n")
    } else {
      cat("  priorfacload:", paste(dim(object$priorfacload), collapse = " x "), "matrix\n")
    }
  }
  if ("restrict" %in% object$supplied_args) {
    cat("  restrict:", object$restrict, "\n")
  }
  cat("----------------------------\n")
  cat("Other:\n")
  cat("  verbose:", object$verbose, "\n")
}



#' Printing method for class mfbvar
#'
#' Method for printing \code{mfbvar} objects.
#'
#' @param x object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @return No return value, called for side effects.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' mod_minn
print.mfbvar <- function(x, ...) {
  ss <- ifelse(x$prior == "ss", "steady-state ", "")
  freq_type <- ifelse(sum(x$freq == "m") == 0, "Quarterly", ifelse(sum(x$freq == "q") == 0, "Monthly", "Mixed-frequency"))
  var_type <- switch(x$variance,
    iw = "Inverse Wishart",
    diffuse = "Diffuse",
    fsv = sprintf("Factor stochastic volatility (%d factors)", x$mfbvar_prior$n_fac),
    csv = "Common stochastic volatility"
  )
  cat(paste0(
    sprintf("%s BVAR with:\n", freq_type), ncol(x$Y), " variables", ifelse(!is.null(x$names_col), paste0(" (", paste(x$names_col, collapse = ", "), ")"), " "),
    "\nPrior: ", x$prior, "\n",
    "\nError covariance matrix: ", var_type, "\n",
    x$n_lags, " lags\n",
    nrow(x$Y), " time periods", ifelse(!is.null(x$names_row), paste0(" (", x$names_row[1], " - ", x$names_row[length(x$names_row)], ")"), " "), "\n", ifelse(is.null(x$n_fcst), "0", x$n_fcst), " periods forecast\n",
    x$n_reps, " draws used in main chain"
  ))
  cat("\n")
}

#' Summary method for class mfbvar
#'
#' Method for summarizing \code{mfbvar} objects.
#'
#' @param object object of class \code{mfbvar}
#' @param ... Currently not in use.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' summary(mod_minn)
summary.mfbvar <- function(object, ...) {
  print(object)
}

#' Plotting methods for posterior mfbvar objects
#'
#' Methods for plotting posterior mfbvar objects.
#' @param x object of class \code{mfbvar_minn} or \code{mfbvar_ss}
#' @param aggregate_fcst Boolean indicating whether forecasts of the latent monthly series should be aggregated to the quarterly frequency.
#' @param plot_start Time period (date or number) to start plotting from. Default is to to use \code{5*n_fcst} time periods if \code{n_fcst} exists, otherwise the entire sample.
#' @param variables Vector of names or positions of variables to include in the plot of variances
#' @param pred_bands Single number (between \code{0.0} and \code{1.0}) giving the coverage level of forecast intervals.
#' @param ss_bands (Steady-state prior only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior steady-state intervals.
#' @param var_bands (\code{varplot} only) Single number (between \code{0.0} and \code{1.0}) giving the coverage level of posterior intervals for the error standard deviations.
#' @param nrow_facet an integer giving the number of rows to use in the facet
#' @param ... Currently not in use.
#' @return A \code{\link[ggplot2]{ggplot}}.
#' @name plot-mfbvar
#' @examples
#' prior_obj <- set_prior(
#'   Y = mf_usa, d = "intercept",
#'   n_lags = 4, n_reps = 20,
#'   n_fcst = 4, n_fac = 1
#' )
#'
#' prior_intervals <- matrix(c(
#'   1, 3,
#'   4, 8,
#'   1, 3
#' ), ncol = 2, byrow = TRUE)
#' psi_moments <- interval_to_moments(prior_intervals)
#' prior_psi_mean <- psi_moments$prior_psi_mean
#' prior_psi_Omega <- psi_moments$prior_psi_Omega
#' prior_obj <- update_prior(prior_obj,
#'   prior_psi_mean = prior_psi_mean,
#'   prior_psi_Omega = prior_psi_Omega
#' )
#'
#' mod_ss <- estimate_mfbvar(prior_obj, prior = "ss", variance = "fsv")
#' plot(mod_ss)
#' varplot(mod_ss)
#' @rdname plot-mfbvar
plot.mfbvar_ss <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ss_bands = 0.95, ...) {
  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    stop("To plot the forecasts, proper dates must be provided in the input data.")
  }
  fcst_start <- lubridate::as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)


  plot_range_names <- fcst_start %m+% months(-nrow(x$Y):(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(nrow(x$Y) - x$n_fcst * 5, 0):nrow(x$Y)
    } else {
      plot_range <- 1:nrow(x$Y)
    }
  } else {
    plot_start <- tryCatch(lubridate::as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):nrow(x$Y)
    } else {
      stop("Unable to convert plot_start to a date.")
    }
  }


  if (is.null(ss_bands)) {
    ss_level <- c(0.025, 0.975)
  } else {
    ss_level <- c(0.5 - ss_bands / 2, 0.5 + ss_bands / 2)
  }
  if (is.null(pred_bands)) {
    pred_level <- c(0.10, 0.90)
  } else {
    pred_level <- c(0.5 - pred_bands / 2, 0.5 + pred_bands / 2)
  }

  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  names_row <- if (is.null(x$names_row)) 1:nrow(x$Y) else x$names_row
  p <- ggplot(mapping = aes(x = time))

  if (x$n_fcst > 0) {
    ss_lower <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
    ss_median <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
    ss_upper <- rbind(x$d[plot_range, , drop = FALSE], x$mfbvar_prior$d_fcst) %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))
  } else {
    ss_lower <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[1]), ncol = x$n_determ))
    ss_median <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = 0.5), ncol = x$n_determ))
    ss_upper <- x$d[plot_range, , drop = FALSE] %*% t(matrix(apply(x$psi, 2, quantile, prob = ss_level[2]), ncol = x$n_determ))
  }

  if (!is.null(x$psi)) {
    ss <- data.frame(expand.grid(time = plot_range[1]:(plot_range[length(plot_range)] + x$n_fcst), variable = names_col),
      lower = c(ss_lower), median = c(ss_median),
      upper = c(ss_upper)
    )
    ss$value <- c(rbind(as.matrix(x$Y[plot_range, ]), matrix(NA, nrow = x$n_fcst, ncol = x$n_vars)))
    ss_excl <- c(rbind(!is.na(as.matrix(x$Y[plot_range, ])), matrix(TRUE, nrow = x$n_fcst, ncol = x$n_vars)))
    # ss <- na.omit(ss)

    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "#bdbdbd"), alpha = 1) +
      geom_line(data = na.omit(ss), aes(y = value))
  }
  if (x$n_fcst > 0) {
    preds <- predict(x, aggregate_fcst = aggregate_fcst, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(
        data.frame(
          variable = names_col[i],
          time = (1:nrow(x$Y))[last_pos[i]],
          fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
          lower = x$Y[last_pos[i], i],
          median = x$Y[last_pos[i], i],
          upper = x$Y[last_pos[i], i]
        ),
        fcst
      )
    }
    fcst <- mutate(fcst, variable = factor(variable, levels = names_col, labels = names_col))
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <- p + geom_ribbon(
      data = fcst, aes(
        ymin = lower, ymax = upper,
        fill = "grey90"
      ), linetype = "dotted", color = "black",
      alpha = 0.75
    ) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = ifelse(aggregate_fcst,
        "Note: The forecasts for the quarterly variables have been aggregated to the quarterly frequency.",
        "Note: The forecasts are for the underlying variable."
      ))
    names_row <- c(names_row, rownames(x$Z_fcst)[-(1:x$n_lags)])
  }

  if (x$n_fcst > 0) {
    p <- p + scale_fill_manual(
      values = c("grey90" = "grey90", "#bdbdbd" = "#bdbdbd"),
      label = c(
        "grey90" = paste0("Prediction (", 100 * (pred_level[2] - pred_level[1]), "%)"),
        "#bdbdbd" = paste0("Steady state (", 100 * (ss_level[2] - ss_level[1]), "%)")
      )
    ) +
      labs(
        fill = "Intervals",
        title = "Forecasts and posterior steady-state intervals",
        y = "Value",
        x = "Time"
      ) +
      guides(fill = guide_legend(override.aes = list(
        fill = c("#bdbdbd", "grey90"),
        linetype = c("blank", "dotted")
      )))
  } else {
    p <- p + scale_fill_manual(
      values = c("#bdbdbd"),
      label = c(paste0("Steady state (", 100 * (ss_level[2] - ss_level[1]), "%)"))
    ) +
      labs(
        fill = "Intervals",
        title = "Steady-state intervals",
        y = "Value",
        x = "Time"
      ) +
      guides(fill = guide_legend(override.aes = list(
        fill = c("#bdbdbd"),
        linetype = c("blank")
      )))
  }


  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position = "bottom")
  breaks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  breaks <- na.omit(breaks)
  if (any(as.numeric(breaks) > plot_range[length(plot_range)])) {
    break_labels <- c(
      as.character(plot_range_names[as.numeric(breaks)[as.numeric(breaks) <= plot_range[length(plot_range)]]]),
      as.character(preds$fcst_date[min(which(preds$time == breaks[as.numeric(breaks) > plot_range[length(plot_range)]]))])
    )
  } else {
    break_labels <- plot_range_names[as.numeric(breaks)]
  }
  p + scale_x_continuous(
    breaks = as.numeric(breaks),
    labels = break_labels
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' @rdname plot-mfbvar
plot.mfbvar_ssng <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                             pred_bands = 0.8, nrow_facet = NULL, ss_bands = 0.95, ...) {
  plot.mfbvar_ss(x,
    aggregate_fcst = aggregate_fcst, plot_start = plot_start,
    pred_bands = pred_bands, nrow_facet = nrow_facet, ss_bands = ss_bands, ...
  )
}

#' @rdname plot-mfbvar
plot.mfbvar_minn <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                             pred_bands = 0.8, nrow_facet = NULL, ...) {
  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    stop("To plot the forecasts, proper dates must be provided in the input data.")
  }
  fcst_start <- lubridate::as_date(rownames(x$Y)[nrow(x$Y)]) %m+% months(1)


  plot_range_names <- fcst_start %m+% months(-nrow(x$Y):(-1))

  lower <- upper <- value <- NULL
  if (is.null(plot_start)) {
    if (x$n_fcst > 0) {
      plot_range <- max(nrow(x$Y) - x$n_fcst * 5, 0):nrow(x$Y)
    } else {
      plot_range <- 1:nrow(x$Y)
    }
  } else {
    plot_start <- tryCatch(lubridate::as_date(plot_start), error = function(cond) cond)
    if (!inherits(plot_start, "error")) {
      if (!(plot_start %in% plot_range_names)) {
        stop(sprintf("The start date, %s, does not match rownames in the data matrix Y.", plot_start))
      }
      plot_range <- (which(plot_range_names == plot_start)):nrow(x$Y)
    } else {
      stop("Unable to convert plot_start to a date.")
    }
  }




  if (is.null(pred_bands)) {
    pred_level <- c(0.10, 0.90)
  } else {
    pred_level <- c(0.5 - pred_bands / 2, 0.5 + pred_bands / 2)
  }

  names_col <- if (is.null(x$names_col)) paste0("x", 1:x$n_vars) else x$names_col
  p <- ggplot(mapping = aes(x = time))


  ss <- data.frame(expand.grid(time = plot_range[1]:(plot_range[length(plot_range)] + x$n_fcst), variable = names_col))
  ss$value <- c(rbind(as.matrix(x$Y[plot_range, ]), matrix(NA, nrow = x$n_fcst, ncol = x$n_vars)))
  ss_excl <- c(rbind(!is.na(as.matrix(x$Y[plot_range, ])), matrix(TRUE, nrow = x$n_fcst, ncol = x$n_vars)))
  # ss <- na.omit(ss)

  p <- p +
    geom_line(data = na.omit(ss), aes(y = value))

  if (x$n_fcst > 0) {
    preds <- predict(x, aggregate_fcst = aggregate_fcst, pred_bands = pred_bands)
    fcst <- preds
    last_pos <- apply(x$Y, 2, function(yy) max(which(!is.na(yy))))
    for (i in seq_along(last_pos)) {
      fcst <- rbind(
        fcst, data.frame(
          variable = names_col[i],
          time = (1:nrow(x$Y))[last_pos[i]],
          fcst_date = preds$fcst_date[1] %m-% months(preds$time[1] - last_pos[i]),
          lower = x$Y[last_pos[i], i],
          median = x$Y[last_pos[i], i],
          upper = x$Y[last_pos[i], i]
        ),
        fcst
      )
    }
    fcst <- mutate(fcst, variable = factor(variable, levels = names_col, labels = names_col))
    fcst <- fcst[!duplicated(fcst[, 1:2]), ]
    p <- p + geom_ribbon(
      data = fcst, aes(
        ymin = lower, ymax = upper,
        fill = "grey90"
      ), linetype = "dotted", color = "black",
      alpha = 0.75
    ) +
      geom_line(data = fcst, aes(y = median), linetype = "dashed") +
      labs(caption = ifelse(aggregate_fcst,
        "Note: The forecasts for the quarterly variables have been aggregated to the quarterly frequency.",
        "Note: The forecasts are for the underlying variable."
      ))
  }

  if (x$n_fcst > 0) {
    p <- p + scale_fill_manual(
      values = c("grey90" = "grey90"),
      label = c("grey90" = paste0("Prediction (", 100 * (pred_level[2] - pred_level[1]), "%)"))
    ) +
      labs(
        fill = "Intervals",
        title = "Forecast intervals",
        y = "Value",
        x = "Time"
      ) +
      guides(fill = guide_legend(override.aes = list(
        fill = c("#bdbdbd"),
        linetype = c("blank")
      )))
  } else {
    p <- p +
      labs(
        y = "Value",
        x = "Time"
      )
  }


  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position = "bottom")
  breaks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks
  breaks <- na.omit(breaks)
  if (any(as.numeric(breaks) > plot_range[length(plot_range)])) {
    break_labels <- c(
      as.character(plot_range_names[as.numeric(breaks)[as.numeric(breaks) <= plot_range[length(plot_range)]]]),
      as.character(preds$fcst_date[min(which(preds$time == breaks[as.numeric(breaks) > plot_range[length(plot_range)]]))])
    )
  } else {
    break_labels <- plot_range_names[as.numeric(breaks)]
  }
  p + scale_x_continuous(
    breaks = as.numeric(breaks),
    labels = break_labels
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot.mfbvar_dl <- function(x, aggregate_fcst = TRUE, plot_start = NULL,
                           pred_bands = 0.8, nrow_facet = NULL, ...) {
  plot.mfbvar_minn(x,
    aggregate_fcst = aggregate_fcst, plot_start = plot_start,
    pred_bands = pred_bands, nrow_facet = nrow_facet, ...
  )
}



#' Plot method for class \code{mfbvar_prior}
#'
#' Method for plotting \code{mfbvar_prior} objects.
#'
#' @param x object of class \code{mfbvar_prior}
#' @param nrow_facet number of rows in facet
#' @param ... Currently not in use.
#' @details The function plots the data. If the prior moments for the steady-state parameters are available in \code{x}, these are included.
#' @return A \code{\link[ggplot2]{ggplot}}.
#' @examples
#' prior_obj <- set_prior(Y = mf_usa, n_lags = 4, n_reps = 20, n_fcst = 4)
#' plot(prior_obj)
plot.mfbvar_prior <- function(x, nrow_facet = NULL, ...) {
  ss_level <- c(0.025, 0.975)

  names_col <- if (is.null(colnames(x$Y))) paste0("x", 1:x$n_vars) else colnames(x$Y)

  if (!is.null(x$d) & !is.null(x$prior_psi_mean) & !is.null(x$prior_psi_Omega)) {
    ss_flag <- TRUE
  } else {
    ss_flag <- FALSE
  }

  if (!is.null(x$d) & !is.null(x$prior_psi_mean)) {
    ssng_flag <- TRUE
  } else {
    ssng_flag <- FALSE
  }

  if (ss_flag) {
    n_determ <- ncol(x$d)
    ss_lower <- x$d %*% t(matrix(qnorm(ss_level[1], x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))
    ss_median <- x$d %*% t(matrix(qnorm(0.5, x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))
    ss_upper <- x$d %*% t(matrix(qnorm(ss_level[2], x$prior_psi_mean, sqrt(diag(x$prior_psi_Omega))), ncol = n_determ))

    ss <- data.frame(expand.grid(time = as.Date(rownames(x$Y)), variable = names_col),
      lower = c(ss_lower), median = c(ss_median),
      upper = c(ss_upper)
    )
  }

  row_names <- tryCatch(as.Date(rownames(x$Y)), error = function(cond) cond)
  if (inherits(row_names, "error")) {
    row_names <- 1:nrow(x$Y)
  }
  plot_df <- data.frame(expand.grid(time = row_names, variable = names_col))
  plot_df$value <- c(as.matrix(x$Y))
  plot_df <- na.omit(plot_df)

  p <- ggplot(mapping = aes(x = time))

  if (ss_flag) {
    p <- p +
      geom_ribbon(data = ss, aes(ymin = lower, ymax = upper, fill = "#bdbdbd"), alpha = 1) +
      scale_fill_manual(
        values = c("#bdbdbd"),
        label = c(paste0("Steady state (", 100 * (ss_level[2] - ss_level[1]), "%)"))
      ) +
      labs(
        fill = "Intervals",
        title = "Prior steady-state intervals"
      ) +
      guides(fill = guide_legend(override.aes = list(
        fill = c("#bdbdbd"),
        linetype = c("blank")
      )))
  }

  p <- p +
    geom_line(data = plot_df, aes(y = value), alpha = 0.75) +
    labs(
      y = "Value",
      x = "Time"
    )

  if (!ss_flag) {
    p <- p +
      labs(title = "Data")
  }

  if (is.null(nrow_facet)) {
    p <- p + facet_wrap(~variable, scales = "free_y")
  } else {
    p <- p + facet_wrap(~variable, scales = "free_y", nrow = nrow_facet)
  }

  p <- p +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}
