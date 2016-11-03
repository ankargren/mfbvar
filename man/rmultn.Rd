\name{rmultn}
\alias{rmultn}
\title{Generate random numbers from a multivariate normal distribution}
\description{
Generates random numbers from the multivariate normal distribution with parameters (m, Sigma).
}
\usage{
rmultn(m, Sigma)
}
\arguments{
\item{m}{A vector of size \code{p}.}
\item{Sigma}{A covariance matrix of size \code{p * p}.}
}
\details{Implemented in C++ using Bartlett decomposition.}
\value{
\item{x}{A vector of size \code{p}.}
}

