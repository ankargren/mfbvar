\name{build_Z}
\alias{build_Z}
\title{Build the \eqn{Z} matrix}
\description{
Builds the \eqn{Z} matrix, which consists of lags of \eqn{z}.
}
\usage{
build_Z(z, n_lags)
}
\arguments{
\item{z}{A matrix of size \code{(n_T + n_lags) * n_vars} of data.}
\item{n_lags}{Number of lags.}
}
\details{The construction is such that \eqn{Z_{1:T}=(z_{1:T}, z_{0:T-1}, ..., z_{2-k:T-k+1})}.}
\value{
\item{Z}{A matrix of size \code{n_T * (n_vars*n_lags)}.}
}

