\name{gibbs_sampler}
\alias{gibbs_sampler}
\title{Gibbs sampler for the MF-BVAR with a steady-state prior}
\description{
A Gibbs sampler for sampling from the posterior distribution of the MF-SS-BVAR.
}
\usage{
gibbs_sampler(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, M_Lambda, n_reps = 10000, n_burnin = 10000)
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
\item{M_Lambda}{A list of length \code{n_T}, where \code{M_Lambda[[t]]} corresponds to \eqn{M_t\Lambda}. The column dimension of each element should be \code{n_lags*n_vars}, but the row dimension may vary.}
\item{n_reps}{The number of post-burn-in replications.}
\item{n_burnin}{The number of draws to discard as burn-in.}
}
\details{The prior covariance of \eqn{\Pi} given \eqn{\Sigma} is \eqn{\Sigma \otimes \Omega_\Pi}.}

