\name{build_DD}
\alias{build_DD}
\title{Build the \eqn{D} matrix}
\description{
Builds the \eqn{D} matrix.
}
\usage{
build_DD(d, n_lags)
}
\arguments{
\item{d}{Matrix of size \code{(n_T + n_lags) * n_determ}.}
\item{n_lags}{Number of lags.}
}
\value{
\item{DD}{A matrix of size \code{n_T * ((n_lags + 1)*n_determ)} where row \code{t} is \eqn{(d_t', -d_{t-1}', \dots, -d_{t-k}')}.}
}

