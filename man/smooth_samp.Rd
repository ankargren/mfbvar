\name{smooth_samp}
\alias{smooth_samp}
\alias{smooth_samp_u3}
\alias{smoothing}

\title{Smooth and sample from the smoothed distribution}
\description{
Smooths and samples from the smoothed distribution \eqn{p(Z_{1:T}|Y_{1:T}, \Theta)}.
}
\usage{
smooth_samp(mZ, mX, lH, lH0 = NULL, mF, mB, mQ, iT, ip, iq, is, h0, P0 = NULL, X0)
smooth_samp_u3(mZ, mX, lH, lH0 = NULL, mF, mB, mQ, iT, ip, iq, is, h0, P0 = NULL, X0)
smoothing(mZ, mX, lH, mF, mB, mQ, iT, ip, iq, is, h0, P0, X0)
}
\arguments{
\item{mZ}{\code{T * n_vars} matrix with the observations (\code{NA} represents missingness)}
\item{mX}{\code{T * n_determ} matrix with the deterministic terms}
\item{lH}{A \code{T}-long list with selection matrices}
\item{lH0}{?}
\item{mF}{\code{(n_vars*n_lags) * (n_vars*n_lags)} matrix containing parameters (companion form)}
\item{mB}{\code{(n_vars*n_lags) * n_determ} matrix containing parameters for the deterministics}
\item{mQ}{\code{(n_vars*n_lags) * (n_vars*n_lags)} matrix whose \code{n_vars*n_vars} top-left corner is the Cholesky decomposition of the error covariance matrix}
\item{iT}{The sample size (sometimes called \code{n_T})}
\item{ip}{The number of variables (\code{n_vars})}
\item{iq}{The companion-form dimension (\code{n_vars*n_lags})}
\item{is}{The number of deterministic terms (\code{n_determ})}
\item{h0}{The initial state (\code{(n_vars*n_lags)*1})}
\item{P0}{The covariance matrix of the initial state (\code{(n_vars*n_lags)*(n_vars*n_lags)})}
\item{X0}{The initial value of the deterministic vector (\code{n_determ * 1})}
}
\value{
For \code{smooth_samp} or \code{smooth_samp_u3}:
\item{mh}{The draws from the posterior distribution.}
\item{mZ}{The draws with the observational structure as in the original data.}
For \code{smoothing}:
\item{mu}{The smoothed errors.}
}
