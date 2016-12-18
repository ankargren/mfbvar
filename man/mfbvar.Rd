\name{mfbvar}
\alias{mfbvar}
\title{MF-BVAR with a steady-state prior}
\description{
The main function for running the MF-SS-BVAR.
}
\usage{
mfbvar(Y, d, d_fcst, Lambda, n_lags, n_burnin, n_reps, prior_mean_pi,
       lambda1, lambda2, prior_nu = NULL, prior_mean_psi, prior_var_psi, ...)
}
\arguments{
\item{Y}{The data matrix of size \code{(n_T + n_lags) * n_vars} with \code{NA} representing missingness. Must correspond well to the specification of \code{M_Lambda}.}
\item{d}{The matrix of size \code{(n_T + n_lags) * n_determ} of deterministic terms.}
\item{d_fcst}{The deterministic terms for the forecasting period.}
\item{Lambda}{The Lambda matrix.}
\item{n_lags}{The number of lags.}
\item{n_burnin}{The number of burn-in replications.}
\item{n_reps}{The number of post-burn-in replications.}
\item{prior_mean_pi}{Matrix of size \code{n_vars * (n_vars*n_lags)} containing the prior for the mean of the dynamic coefficients.}
\item{lambda1}{The overall tightness.}
\item{lambda2}{The lag decay.}
\item{prior_nu}{The prior degrees of freedom in the inverse Wishart distribution for Sigma. If \code{NULL}, it is set to \code{n_vars+2} (the lowest value to yield a proper prior).}
\item{prior_mean_psi}{Vector of length \code{n_determ*n_vars} with the prior means of the steady-state parameters.}
\item{prior_var_psi}{Matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)} with the prior covariance of the steady-state parameters.}
\item{...}{Other calls to the main run (post-burn-in) of \code{\link{gibbs_sampler}}.}
}
\seealso{
\code{\link{gibbs_sampler}}
}
