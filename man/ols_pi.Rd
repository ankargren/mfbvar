\name{ols_pi}
\alias{ols_pi}
\alias{ols_s}
\title{OLS functions}
\description{
Compute OLS estimates and error covariance.
}
\usage{
ols_pi(X, Y)
ols_s(X, Y, Pi)
}
\arguments{
\item{X}{The regressor matrix.}
\item{Y}{The matrix of dependent variables.}
\item{Pi}{The coefficient matrix (as returned by \code{ols_pi}).}
}
\value{
One of:
\item{pi_sample}{Estimated coefficients.}
\item{s_sample}{The sum of squared residuals matrix.}
}

