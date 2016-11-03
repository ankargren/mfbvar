\name{posterior_psi}
\alias{posterior_psi}
\alias{posterior_psi_omega}
\title{Compute posterior moments of the steady-state parameters}
\description{
Computes the mean and variance of the conditional posterior distribution of the steady-state parameters.
}
\usage{
posterior_psi_omega(U, D, sigma, prior_psi_omega)
posterior_psi(U, D, sigma, prior_psi_omega, psi_omega, Y_tilde, prior_psi)
}
\arguments{
\item{U}{The \eqn{U} matrix, of size \code{(n_vars*n_determ*(n_lags+1)) * (n_vars*n_determ)}. This can be obtained using \code{\link{build_U}}.}
\item{D}{The \eqn{D} matrix of size \code{n_T * (n_determ*(n_lags+1))} with the deterministic terms and lags thereof. Can be obtained using \code{\link{build_DD}}.}
\item{sigma}{The error covariance matrix of size \code{n_vars * n_vars}.}
\item{prior_psi_omega}{The prior for \eqn{\Omega_\psi}, a matrix of size \code{(n_determ*n_vars) * (n_determ*n_vars)}.}
\item{psi_omega}{The posterior of \eqn{\Omega_\psi} (i.e. obtained using \code{\link{posterior_psi_omega}}).}
\item{Y_tilde}{The lag-corrected data matrix (with no missing values) of size \code{n_T * n_vars}.}
\item{prior_psi}{The prior mean of \eqn{\psi}, a vector (or matrix) of length \code{n_determ * n_vars}.}
}
\details{}
\value{The functions return one of:
\item{psi}{The posterior mean (from \code{\link{posterior_psi}})}
\item{psi_omega}{The posterior variance (from \code{\link{posterior_psi_omega}})}
}

