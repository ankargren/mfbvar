\name{mdd}
\alias{mdd}
\alias{mdd1}
\title{Marginal data density}
\description{
  Estimate the marginal data density.
}
\usage{
  mdd(mfbvar_obj, p_trunc)
  mdd1(mfbvar_obj)
}
\arguments{
  \item{mfbvar_obj}{The return from \code{\link{mfbvar}}.}
  \item{p_trunc}{\code{1-p_trunc} is the degree of truncation.}
}
\value{A list with components being \code{n_reps}-long vectors.
  \item{pi_sigma_posterior}{Posterior of Pi and Sigma.}
  \item{data_likelihood}{The likelihood.}
  \item{pi_sigma_prior}{Prior of Pi and Sigma.}
  \item{psi_prior}{Prior of psi.}
  \item{psi_truncated}{The truncated psi pdf.}
}
\details{
See the example for computation of the marginal data density.
}
\examples{
\dontrun{
mdd_out <- mdd(mfbvar_obj, 0.9)
mdd_num <- mdd_out$pi_sigma_posterior * mdd_out$psi_truncated
mdd_den <- mdd_out$data_likelihood * mdd_out$pi_sigma_prior * mdd_out$psi_prior
mdd_res <- 1/sum(mdd_num/mdd_den)
}
}

