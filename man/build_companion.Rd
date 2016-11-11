\name{build_companion}
\alias{build_companion}
\title{Build the companion matrix for the dynamic parameters}
\description{
Builds the parameter matrix of dynamic coefficients for the companion form representation.
}
\usage{
build_companion(Pi, n_vars, n_lags)
}
\arguments{
\item{Pi}{The parameter matrix of size \code{n_vars * (n_vars*n_lags)}.}
\item{n_vars}{Number of variables.}
\item{n_lags}{Number of lags.}
}
\value{
\item{Pi_comp}{The companion form matrix of size \code{(n_vars*n_lags) * (n_vars*n_lags)}.}
}

