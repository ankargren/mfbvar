\name{dfuns}
\alias{dmultn}
\alias{dnorm_trunc}
\alias{dnorminvwish}
\title{Density functions}
\description{
  Functions for evaluations of pdfs.
}
\usage{
  dmultn(x, m, Sigma)
  dnorm_trunc(x, m, V_inv, d, p_trunc, chisq_val)
  dnorminvwish(X, Sigma, M, P, S, v)
}
\arguments{
  \item{x}{A vector of size \code{p}.}
  \item{X}{Matrix of size \code{p * q}.}
  \item{m}{The mean vector of size \code{p}.}
  \item{M}{The mean matrix of size \code{p * q}.}
  \item{Sigma}{The covariance matrix.}
  \item{V_inv}{The inverse of the covariance matrix (for \code{dnorm_trunc}) of size \code{d * d}.}
  \item{d}{The length of the vector (for \code{dnorm_trunc}).}
  \item{p_trunc}{\code{1-p_trunc} is the degree of truncation (i.e. \code{p_trunc=1} is no truncation).}
  \item{chisq_val}{The value in the corresponding chi-square distribution; if the normal quadratic form exceeds this, the pdf is 0.}
  \item{P}{\code{p * p} covariance matrix.}
  \item{S}{\code{q * q} scale matrix for the inverse Wishart.}
  \item{v}{The degrees of freedom for the inverse Wishart.}
}
\value{
  The value of the pdf.
}

