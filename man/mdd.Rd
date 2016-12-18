\name{mdd}
\alias{mdd}
\title{Marginal data density}
\description{
  Estimate the marginal data density.
}
\usage{
  mdd(mfbvar_obj, p_trunc)
}
\arguments{
  \item{mfbvar_obj}{The return from \code{\link{mfbvar}}.}
  \item{p_trunc}{\code{1-p_trunc} is the degree of truncation.}
}
\value{
  A vector, the mean of which is the estimated marginal data density.
}

