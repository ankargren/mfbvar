\name{build_Y_tilde}
\alias{build_Y_tilde}
\title{Build the lag-corrected data matrix}
\description{
Builds the \eqn{\tilde{Y}=\Pi(L)Y} matrix.
}
\usage{
build_Y_tilde(Pi, z)
}
\arguments{
\item{Pi}{A matrix of size \code{n_vars * (n_vars*n_lags)} of dynamic coefficients.}
\item{z}{A matrix of size \code{(n_T + n_lags) * n_vars}.}
}
\details{Note that \code{z} does not contain missing values; at this point, the missing values have been replaced by values drawn using the simulation smoother.}
\value{
\item{Y_tilde}{A matrix of size \code{n_t * n_vars}.}
}

