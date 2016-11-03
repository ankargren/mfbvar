\name{rmatn}
\alias{rmatn}
\title{Generate random numbers from a matricvariate normal distribution}
\description{
Generates random numbers from the matricvariate normal distribution with parameters (M, Q, P).
}
\usage{
rmatn(M, Q, P)
}
\arguments{
\item{M}{A matrix of size \code{p * q} with the means.}
\item{Q}{A matrix of size \code{q * q} with the first covariance matrix.}
\item{P}{A matrix of size \code{p * p} with the first covariance matrix.}
}
\details{Implemented in C++.}
\value{
\item{X}{A matrix of size \code{p * q}.}
}

