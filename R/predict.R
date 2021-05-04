#' Predict method for class \code{mfbvar}
#'
#' Method for predicting \code{mfbvar} objects.
#'
#' @param object object of class mfbvar
#' @param aggregate_fcst If forecasts of quarterly variables should be aggregated back to the quarterly frequency.
#' @param pred_bands The level of the probability bands for the forecasts.
#' @param ... Currently not in use.
#' @details Note that this requires that forecasts were made in the original \code{mfbvar} call.
#' @return A \code{\link[tibble]{tibble}} with columns:
#' \describe{\item{\code{variable}}{Name of variable}
#' \item{\code{time}}{Time index}
#' \item{\code{fcst_date}}{Date of forecast}}
#' If the argument \code{pred_bands} is given as a numeric value between 0 and 1, the returned tibble also includes columns:
#' \describe{\item{\code{lower}}{The \code{(1-pred_bands)/2} lower quantiles of the predictive distributions}
#' \item{\code{median}}{The medians of the predictive distributions}
#' \item{\code{upper}}{The \code{(1+pred_bands)/2} upper quantiles of the predictive distributions}}
#' If \code{pred_bands} \code{NULL} or \code{NA}, the returned tibble also includes the columns:
#' \describe{\item{\code{fcst}}{MCMC samples from the predictive distributions}
#' \item{\code{iter}}{Iteration indexes for the MCMC samples}}
#' @examples
#' prior_obj <- set_init(Y = mf_usa, n_lags = 4, n_burnin = 20, n_reps = 20,
#' n_fcst = 4)
#' prior_obj <- set_prior_minn(prior_obj)
#' mod_minn <- estimate_mfbvar(prior_obj, prior = "minn")
#' predict(mod_minn)
predict.mfbvar <- function(object, aggregate_fcst = TRUE, pred_bands = 0.8, ...) {
  end_month <- FALSE
  if (object$n_fcst == 0) {
    stop("No forecasts exist in the provided object.")
  }

  final_h <- c(unlist(apply(object$Y, 2, function(x) Position(is.na, x, nomatch = nrow(object$Y) + 1)))[object$mfbvar_prior$freq == object$mfbvar_prior$freqs[2]]) - 1
  final_l <- min(apply(object$Y[, object$mfbvar_prior$freq == object$mfbvar_prior$freqs[1], drop = FALSE], 2, function(x) max(which(!is.na(x)))))
  final_non_na <- min(c(
    final_h,
    final_l
  ))
  final_fcst <- object$n_lags - (nrow(object$Y) - final_non_na) + 1
  if (final_fcst >= 1) {
    incl_fcst <- final_fcst:(object$n_lags + object$n_fcst)
  } else {
    incl_fcst <- 1:(object$n_lags + object$n_fcst)
  }

  n_h <- sum(object$mfbvar_prior$freq == object$mfbvar_prior$freqs[2])
  n_l <- sum(object$mfbvar_prior$freq == object$mfbvar_prior$freqs[1])
  n_vars <- n_h + n_l

  tmp <- tryCatch(lubridate::ymd(rownames(object$Y)[nrow(object$Y)]), warning = function(cond) cond)

  if (!("w" %in% object$freq) && !inherits(tmp, "warning")) {
    final_est <- lubridate::ymd(rownames(object$Y)[nrow(object$Y)])
    fcst_start <- final_est %m+% months(1)
    if (lubridate::days_in_month(final_est) == lubridate::day(final_est)) {
      end_month <- TRUE
    }

    ret_names <- fcst_start %m+% months((-(length(incl_fcst) - object$n_fcst)):(object$n_fcst - 1))
    ret_names_q <- fcst_start %m+% months((-(object$n_lags)):(object$n_fcst - 1))

    if (end_month) {
      ret_names <- lubridate::ceiling_date(ret_names, unit = "months") - lubridate::days(1)
      ret_names_q <- lubridate::ceiling_date(ret_names_q, unit = "months") - lubridate::days(1)
    }
  } else {
    final_est <- nrow(object$Y)
    fcst_start <- final_est + 1

    ret_names <- fcst_start + (-(length(incl_fcst) - object$n_fcst)):(object$n_fcst - 1)


    if (aggregate_fcst) {
      if ("w" %in% object$freq) {
        warning("Because of ambiguities, forecasts are not aggregated when weekly data are included.")
      } else {
        stop("Dates must be provided to aggregate latent monthly forecasts to the quarterly frequency.")
      }
    }
  }



  if (!aggregate_fcst) {
    fcst_collapsed <- tibble(
      variable = rep(rep(as.character(object$names_col), each = length(incl_fcst)), object$n_reps / object$n_thin),
      iter = rep(1:(object$n_reps / object$n_thin), each = n_vars * length(incl_fcst)),
      fcst = c(object$Z_fcst[incl_fcst, , ]),
      fcst_date = rep(ret_names, n_vars * object$n_reps / object$n_thin),
      freq = rep(rep(object$freq, each = length(incl_fcst)), object$n_reps / object$n_thin),
      time = rep(nrow(object$Y) + object$n_fcst - max(incl_fcst) + incl_fcst, n_vars * object$n_reps / object$n_thin)
    )
  } else {
    fcst_collapsed <- tibble(
      variable = rep(rep(object$names_col[1:n_h], each = length(incl_fcst)), object$n_reps / object$n_thin),
      iter = rep(1:(object$n_reps / object$n_thin), each = n_h * length(incl_fcst)),
      fcst = c(object$Z_fcst[incl_fcst, 1:n_h, ]),
      fcst_date = rep(as.Date(as.character(ret_names)), n_h * object$n_reps / object$n_thin),
      freq = rep(rep(rep(object$mfbvar_prior$freqs[2], n_h), each = length(incl_fcst)), object$n_reps / object$n_thin),
      time = rep(nrow(object$Y) + object$n_fcst - max(incl_fcst) + incl_fcst, n_h * object$n_reps / object$n_thin)
    ) %>%
      transmute(
        variable = variable,
        iter = iter,
        year = year(fcst_date),
        quarter = quarter(fcst_date),
        fcst_date = fcst_date,
        fcst = fcst,
        freq = freq,
        time = time
      )
    if (aggregate_fcst) {
      n_Lambda <- ncol(object$Lambda_) / nrow(object$Lambda_)
      fcst_agg_required <- final_l + 3 - n_Lambda + 1
      fcst_included <- nrow(object$Y) - object$n_lags + 1
      fcst_agg_missing <- max(c(fcst_included - fcst_agg_required, 0))
      fcst_q <- array(0, dim = c(dim(object$Z_fcst)[1] + max(c(fcst_agg_missing, 0)), n_l, object$n_reps / object$n_thin))
      if (fcst_agg_required < fcst_included) {
        ret_names_q <- c(
          ret_names_q[1] %m+% months((-fcst_agg_missing):(-1)),
          ret_names_q
        )
        fcst_q[1:fcst_agg_missing, , ] <- object$Z[fcst_agg_required:(fcst_included - 1), object$mfbvar_prior$freq == object$mfbvar_prior$freqs[1], , drop = FALSE]
      } else {
        if (nrow(fcst_q) > length(ret_names_q)) {
          ret_names_q <- c(
            ret_names_q[1] %m+% months((-(nrow(fcst_q) - length(ret_names_q))):(-1)),
            ret_names_q
          )
        }
      }
      fcst_q[(fcst_agg_missing + 1):nrow(fcst_q), , ] <- object$Z_fcst[, object$mfbvar_prior$freq == object$mfbvar_prior$freqs[1], , drop = FALSE]
      rownames(fcst_q) <- as.character(ret_names_q)

      end_of_quarter <- which(lubridate::month(ret_names_q) %% 3 == 0)
      end_of_quarter <- end_of_quarter[end_of_quarter >= n_Lambda]
      agg_fun <- function(fcst_q, Lambda_, end_of_quarter) {
        fcst_q_agg <- array(0, dim = c(length(end_of_quarter), dim(fcst_q)[2:3]))
        for (i in 1:(object$n_reps / object$n_thin)) {
          Z_i <- matrix(fcst_q[, , i], nrow = nrow(fcst_q), ncol = ncol(fcst_q))
          for (j in 1:length(end_of_quarter)) {
            Z_ij <- matrix(t(Z_i[(((-n_Lambda + 1):0) + end_of_quarter[j]), , drop = FALSE]), ncol = 1)
            fcst_q_agg[j, , i] <- Lambda_ %*% Z_ij
          }
        }
        return(fcst_q_agg)
      }

      fcst_q_agg <- agg_fun(fcst_q, object$Lambda_, end_of_quarter)

      fcst_quarterly <- tibble(
        variable = rep(rep(object$names_col[(n_h + 1):n_vars], each = nrow(fcst_q_agg)), object$n_reps / object$n_thin),
        iter = rep(1:(object$n_reps / object$n_thin), each = n_l * nrow(fcst_q_agg)),
        fcst = c(fcst_q_agg),
        fcst_date = rep(ret_names_q[end_of_quarter], n_l * object$n_reps / object$n_thin),
        freq = rep(rep(rep(object$mfbvar_prior$freqs[1], n_l), each = nrow(fcst_q_agg)), object$n_reps / object$n_thin),
        time = rep(seq(final_l + 3, by = 3, length.out = nrow(fcst_q_agg)), n_l * object$n_reps / object$n_thin)
      ) %>%
        transmute(
          variable = variable,
          iter = iter,
          year = year(fcst_date),
          quarter = quarter(fcst_date),
          fcst_date = fcst_date,
          fcst = fcst,
          freq = freq,
          time = time
        )
    } else {
      fcst_quarterly <- tibble(
        variable = rep(rep(object$names_col[(n_h + 1):n_vars], each = length(incl_fcst)), object$n_reps / object$n_thin),
        iter = rep(1:(object$n_reps / object$n_thin), each = n_l * length(incl_fcst)),
        fcst = c(object$Z_fcst[incl_fcst, (n_h + 1):n_vars, ]),
        fcst_date = rep(ret_names, n_l * object$n_reps / object$n_thin),
        freq = rep(rep(rep(object$mfbvar_prior$freqs[1], n_l), each = length(incl_fcst)), object$n_reps / object$n_thin),
        time = rep(nrow(object$Y) + object$n_fcst - max(incl_fcst) + incl_fcst, n_l * object$n_reps / object$n_thin)
      ) %>%
        transmute(
          variable = variable,
          iter = iter,
          year = year(fcst_date),
          quarter = quarter(fcst_date),
          fcst_date = fcst_date,
          fcst = fcst,
          freq = freq,
          time = time
        )
    }

    fcst_collapsed <- bind_rows(fcst_collapsed, fcst_quarterly)
  }

  if (!is.null(pred_bands) && !is.na(pred_bands)) {
    pred_quantiles <- c(0.5 - pred_bands / 2, 0.5, 0.5 + pred_bands / 2)
    fcst_collapsed <- group_by(fcst_collapsed, variable, time, fcst_date) %>%
      summarize(
        lower = quantile(fcst, prob = pred_quantiles[1], names = FALSE),
        median = quantile(fcst, prob = pred_quantiles[2], names = FALSE),
        upper = quantile(fcst, prob = pred_quantiles[3], names = FALSE),
        .groups = "keep"
      ) %>%
      ungroup()
  } else {
    fcst_collapsed <- fcst_collapsed[, c("variable", "time", "fcst_date", "fcst", "iter")]
  }

  return(fcst_collapsed)
}

predict.sfbvar <- function(object, pred_bands = 0.8, ...) {
  end_period <- FALSE
  sf_type <- unique(object$mfbvar_prior$freq)
  if (object$n_fcst == 0) {
    stop("No forecasts exist in the provided object.")
  }
  if (object$n_fcst > 0) {
    tmp <- tryCatch(lubridate::ymd(rownames(object$Y)[nrow(object$Y)]), warning = function(cond) cond)
    if (inherits(tmp, "warning")) {
      stop("To summarize the forecasts, proper dates must be provided in the input data.")
    } else {
      final_est <- lubridate::ymd(rownames(object$Y)[nrow(object$Y)])
      if (sf_type == "m") {
        fcst_start <- final_est %m+% months(1)
      } else {
        fcst_start <- final_est %m+% months(3)
      }
      if (lubridate::days_in_month(final_est) == lubridate::day(final_est)) {
        end_period <- TRUE
      }
    }
  }

  final_fcst <- object$n_lags + 1
  if (final_fcst >= 1) {
    incl_fcst <- final_fcst:(object$n_lags + object$n_fcst)
  } else {
    incl_fcst <- 1:(object$n_lags + object$n_fcst)
  }

  if (sf_type == "m") {
    ret_names <- fcst_start %m+% months((-(length(incl_fcst) - object$n_fcst)):(object$n_fcst - 1))
  } else {
    ret_names <- fcst_start %m+% months(3 * ((-(length(incl_fcst) - object$n_fcst)):(object$n_fcst - 1)))
  }

  if (end_period) {
    ret_names <- lubridate::ceiling_date(ret_names, unit = "months") - lubridate::days(1)
  }

  fcst_collapsed <- tibble(
    variable = rep(rep(object$names_col, each = length(incl_fcst)), object$n_reps / object$n_thin),
    iter = rep(1:(object$n_reps / object$n_thin), each = object$n_vars * length(incl_fcst)),
    fcst = c(object$Z_fcst[incl_fcst, , ]),
    fcst_date = rep(as.Date(as.character(ret_names)), object$n_vars * object$n_reps / object$n_thin),
    freq = rep(rep(object$mfbvar_prior$freq, each = length(incl_fcst)), object$n_reps / object$n_thin),
    time = rep(nrow(object$Y) + object$n_fcst - max(incl_fcst) + incl_fcst, object$n_vars * object$n_reps / object$n_thin)
  ) %>%
    transmute(
      variable = variable,
      iter = iter,
      year = year(fcst_date),
      quarter = quarter(fcst_date),
      fcst_date = fcst_date,
      fcst = fcst,
      freq = freq,
      time = time
    )


  if (!is.null(pred_bands) && !is.na(pred_bands)) {
    pred_quantiles <- c(0.5 - pred_bands / 2, 0.5, 0.5 + pred_bands / 2)
    fcst_collapsed <- group_by(fcst_collapsed, variable, time, fcst_date) %>%
      summarize(
        lower = quantile(fcst, prob = pred_quantiles[1], names = FALSE),
        median = quantile(fcst, prob = pred_quantiles[2], names = FALSE),
        upper = quantile(fcst, prob = pred_quantiles[3], names = FALSE)
      ) %>%
      ungroup()
  }


  return(fcst_collapsed)
}
