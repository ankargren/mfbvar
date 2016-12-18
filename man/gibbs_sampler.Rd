\name{gibbs_sampler}
\alias{gibbs_sampler}
\alias{gibbs_sampler_qf}
\title{Gibbs sampler for the MF-BVAR with a steady-state prior}
\description{
A Gibbs sampler for sampling from the posterior distribution of the MF-SS-BVAR.
}
\usage{
gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, n_reps, n_fcst = NULL, lH, check_roots = TRUE,
                          init_Pi = NULL, init_Sigma = NULL, init_psi = NULL,
                          init_Z = NULL, d_fcst = NULL, smooth_state = FALSE)
}
\arguments{
\item{prior_pi}{Matrix of size \code{n_vars * (n_vars*n_lags)} containing the prior for the mean of the dynamic coefficients.}
\item{prior_pi_omega}{Matrix of size \code{(n_vars*n_lags)* (n_vars*n_lags)} containing the prior for (part of) the prior covariance of the dynamic coefficients.}
\item{prior_nu}{The prior degrees of freedom.}
\item{prior_s}{The prior for \eqn{\Sigma}.}
\item{prior_psi}{Vector of length \code{n_determ*n_vars} with the prior means of the steady-state parameters.}
\item{prior_psi_omega}{Matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)} with the prior covariance of the steady-state parameters.}
\item{Y}{The data matrix of size \code{(n_T + n_lags) * n_vars} with \code{NA} representing missingness. Must correspond well to the specification of \code{M_Lambda}.}
\item{d}{The matrix of size \code{(n_T + n_lags) * n_determ} of deterministic terms.}
\item{n_reps}{The number of post-burn-in replications.}
\item{n_fcst}{The number of periods to forecast.}
\item{lH}{A list of length \code{n_T}, where \code{M_Lambda[[t]]} corresponds to \eqn{M_t\Lambda}. The column dimension of each element should be \code{n_lags*n_vars}, but the row dimension may vary.}
\item{check_roots}{Logical, if roots of the companion matrix are to be checked to ensure stationarity.}
\item{init_Pi}{Matrix to initialize the dynamic coefficients.}
\item{init_Sigma}{Matrix to initialize the error covariance.}
\item{init_psi}{Matrix to initialize the steady-state parameters.}
\item{init_Z}{Matrix to initialize the latent state.}
\item{d_fcst}{The deterministic terms for the forecasting period.}
\item{smooth_state}{Logical, if \code{TRUE} then the smoothed estimates of the latent states are also returned.}
}
\details{The prior covariance of \eqn{\Pi} given \eqn{\Sigma} is \eqn{\Sigma \otimes \Omega_\Pi}.

The function \code{gibbs_sampler_qf} is a temporary version of \code{gibbs_sampler} which is customized for quarterly data (for improved speed). Its purpose is to accomomdate conditional forecasting.}

