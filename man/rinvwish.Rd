\name{rinvwish}
\alias{rinvwish}
\title{Generate random numbers from the inverse Wishart distribution}
\description{
Generates random numbers from the inverse Wishart distribution with parameters (v, S).
}
\usage{
rinvwish(v, S)
}
\arguments{
\item{v}{The degrees of freedom.}
\item{S}{A matrix of size \code{p * p}.}
}
\details{Implemented in C++ using Bartlett decomposition.}
\value{
\item{X}{A matrix of size \code{p * p}.}
}

