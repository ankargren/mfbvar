# The below loop just gets the error variances from AR(4) regressions

compute_error_variances <- function(Y) {
  n_vars <- ncol(Y)
  error_variance <- rep(NA, n_vars)
  for (i in 1:n_vars) {
    success <- NULL
    init_order <- 4
    while(is.null(success)) {
      error_variance[i] <- tryCatch(arima(na.omit(Y[,i]), order = c(init_order, 0, 0), method = "ML")$sigma2,
                                    error = function(cond) NA)
      if (!is.na(error_variance[i])) {
        success <- 1
      } else {
        init_order <- init_order - 1
        if (init_order < 1) {
          error_variance[i] <- var(na.omit(Y[,i]))
        }
      }
    }
  }
  return(error_variance)
}
