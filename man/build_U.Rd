\name{build_U}
\alias{build_U}
\alias{build_U_cpp}
\title{Build the \eqn{U} matrix}
\description{
Builds the \eqn{U} matrix, which consists of \eqn{\Pi}.
}
\usage{
build_U(Pi, n_determ)
build_U_cpp(Pi, n_determ, n_vars, n_lags)
}
\arguments{
\item{Pi}{A matrix of size \code{n_vars * (n_vars*n_lags)} of dynamic coefficients.}
\item{n_determ}{Number of deterministic terms.}
\item{n_vars}{Number of variables.}
\item{n_lags}{Number of lags.}
}
\details{\code{build_U} is implemented in R, and \code{build_U_cpp} in C++. Note that they have a different number of arguments.}
\value{
\item{U}{A matrix of size \code{(n_vars*n_determ*(n_lags+1)) * (n_vars*n_determ)}.}
}

