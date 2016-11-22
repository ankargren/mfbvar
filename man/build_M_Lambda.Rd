\name{build_M_Lambda}
\alias{build_M_Lambda}
\title{Build the \eqn{M_t\Lambda} matrices}
\description{
Builds the \eqn{M_t\Lambda} matrices and stores them in a list.
}
\usage{
build_M_Lambda(Y, Lambda, n_vars, n_lags, n_T)
}
\arguments{
\item{Y}{The data matrix of size \code{n_T * n_vars} containing some missing values (\code{NA}).}
\item{Lambda}{Matrix of size \code{n_vars * (n_vars*n_lags)} which describes the aggregation scheme.}
\item{n_vars}{Number of variables.}
\item{n_lags}{Number of lags.}
\item{n_T}{Number of time periods.}
}
\details{The element \code{M_Lambda[[t]]} corresponds to \eqn{M_t\Lambda}. Currently, if element \code{i} of \code{Y[t, ]} is \code{NA}, then row \code{i} of \code{M_Lambda[[t]]} is all \code{NA}.}
\value{
\item{M_Lambda}{A list of length \code{n_T}.}
}

