\name{build_demeaned_z}
\alias{build_demeaned_z}
\alias{build_demeaned_z_cpp}
\title{Demean the data matrix}
\description{
Builds the \eqn{\tilde{z}_{1:T}} matrix, i.e. the demeaned \eqn{z_{1:T}}.
}
\usage{
build_demeaned_z(z, psi, d)
build_demeaned_z_cpp(z, psi, d, n_T, n_vars)
}
\arguments{
\item{z}{A matrix of size \code{n_T * n_vars} of data.}
\item{psi}{A matrix of size \code{(n_vars*n_determ) * 1} of vectorized steady-state parameters.}
\item{d}{A matrix of size \code{(n_T + n_vars) * n_determ} of deterministic data.}
\item{n_T}{The length of the time series (i.e. \code{nrow(z)}).}
\item{n_vars}{The number of time series (i.e. \code{ncol(z)}).}
}
\details{\code{build_demeaned_z} is implemented in R, and \code{build_demeaned_z_cpp} in C++. Note that they have a different number of arguments.}
\value{
\item{demeaned_z}{A matrix of size \code{n_T * n_vars}.}
}

