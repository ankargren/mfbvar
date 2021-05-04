#' mfbvar: A package for mixed-frequency Bayesian vector autoregressive (VAR) models.
#'
#' The mfbvar package makes estimation of Bayesian VARs with a mix of monthly and quarterly data
#' simple. The prior for the regression parameters is normal with Minnesota-style prior moments.
#' The package supports either an inverse Wishart prior for the error covariance matrix, yielding a
#' standard normal-inverse Wishart prior, or a time-varying error covariance matrix by means of a factor
#' stochastic volatility model through the \code{\link[factorstochvol]{factorstochvol-package}} package.
#'
#' @section Specifying the prior:
#' The prior of the VAR model is specified using the function \code{\link{set_prior}} or a selected choice of its wrapper functions. The function creates a prior object, which can be further updated using the wrappers or itself in an object-oriented manner. The model can be estimated using the steady-state prior, which requires the prior moments of the steady-state parameters. The function \code{\link{interval_to_moments}} is a helper function for obtaining these from prior intervals.
#'
#' @section Estimating the model:
#' The model is estimated using the function \code{\link{estimate_mfbvar}}. The error covariance matrix
#' is given an inverse Wishart or diffuse prior or modeled using factor or common stochastic volatility. If the inverse Wishart prior is used,
#' \code{\link{mdd}} can be used to estimate to the marginal data density (marginal likelihood).
#'
#' @section Processing the output:
#' Plots of the output can be obtained from calling the generic function \code{plot} (see
#' \code{\link{plot-mfbvar}}). If factor stochastic volatility is used, the time-varying
#' standard deviations can be plotted using \code{\link{varplot}}. Predictions can be obtained
#' from \code{\link{predict.mfbvar}}.
#'
#'
#' @docType package
#' @name mfbvar
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", "obj", "prior_type", "lower", "upper", "value",
    "variable", "iter", "fcst_date", "fcst", "freq", "prior_Pi_AR1"
  ))
}
