\name{smoothing}
\alias{smoother}
\alias{simulation_smoother}
\alias{generate_mhh}

\title{Smooth and sample from the smoothed distribution}
\description{
Smooths and samples from the smoothed distribution \eqn{p(Z_{1:T}|Y_{1:T}, \Theta)}.
}
\usage{
smoother(mZ, Lambda, mF, mQ, iT, ip, iq, h0, P0)
simulation_smoother(mZ, Lambda, mF, mQ, iT, ip, iq, h0, P0)
generate_mhh(mZ, Lambda, mF, mQ, iT, ip, iq, h0, P0)
}
\arguments{
\item{mZ}{\code{T * n_vars} matrix with the observations (\code{NA} represents missingness)}
\item{Lambda}{The Lambda selection matrix (size \code{n_vars * (n_vars*n_lags)}).}
\item{mF}{\code{(n_vars*n_lags) * (n_vars*n_lags)} matrix containing parameters (companion form)}
\item{mQ}{\code{(n_vars*n_lags) * (n_vars*n_lags)} matrix whose \code{n_vars*n_vars} top-left block is the Cholesky decomposition of the error covariance matrix}
\item{iT}{The sample size (sometimes called \code{n_T})}
\item{ip}{The number of variables (\code{n_vars})}
\item{iq}{The companion-form dimension (\code{n_vars*n_lags})}
\item{h0}{The initial state (\code{(n_vars*n_lags)*1})}
\item{P0}{The covariance matrix of the initial state (\code{(n_vars*n_lags)*(n_vars*n_lags)})}
}
\value{
For \code{simulation_smoother}:
\item{}{The draw from the posterior distribution.}
For \code{smoother}:
\item{}{The smoothed states.}
For \code{generate_mhh}:
\item{}{Generated (pseudo-)state vector.}
}
\details{
Implemented in C++.
}
