\name{ols_initialization}
\alias{ols_initialization}
\title{Initialize Gibbs sampler using OLS}
\description{
Initializes the Gibbs sampler using OLS.
}
\usage{
ols_initialization(Z, d, n_lags)
}
\arguments{
\item{z}{Matrix of size \code{(n_T + n_lags) * n_vars} of data.}
\item{d}{Matrix of size \code{(n_T + n_lags) * n_determ} of deterministic data}
\item{n_lags}{Number of lags.}
}
\value{
A list with components:
\item{Gam}{A matrix of size \code{n_vars * (n_vars*n_lags +n_determ)} of estimated parameters.}
\item{S}{Estimated error covariance matrix.}
}

