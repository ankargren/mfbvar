#' MCMC sampler
#'
#' \code{mcmc_sampler} is a generic function for deciding which specific MCMC
#' algorithm to dispatch to. It is called internally.
#'
#' @param x argument to dispatch on (of class \code{prior_obj})
#' @param ... additional named arguments passed on to the methods

mcmc_sampler <- function(x, ...) {
  UseMethod("mcmc_sampler")
}
