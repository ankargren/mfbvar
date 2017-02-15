\name{eval}
\alias{eval_Pi_Sigma_RaoBlack}
\alias{eval_psi_MargPost}
\title{Functions for evaluating posteriors}
\description{
The functions are helper functions for evaluating posterior densities.
}
\usage{
eval_Pi_Sigma_RaoBlack(Z, d, post_psi, post_pi_mean, post_Sigma, nu, prior_pi,
                       prior_pi_omega, prior_s, n_vars, n_lags, n_reps)
eval_psi_MargPost(Pi, Sigma, Z, post_psi, prior_psi, prior_psi_omega, D,
                  n_determ, n_vars, n_lags, n_reps)
}
\arguments{
\item{Pi}{The array of dynamic coefficient matrices.}
\item{Sigma}{The array of error covariance matrices.}
\item{post_psi}{The posterior center of psi.}
\item{Z}{The array of latent states.}
\item{post_pi_mean}{Posterior center of dynamic coefficients.}
\item{post_Sigma}{Posterior center of error covariance.}
\item{nu}{Posterior degrees of freedom.}
\item{prior_s}{The prior for Sigma.}
\item{prior_pi}{The prior mean for Pi.}
\item{prior_pi_omega}{The prior covariance matrix for Pi.}
\item{prior_psi}{The prior for the mean of psi.}
\item{prior_psi_omega}{The prior covariance matrix for psi.}
\item{D}{The \code{D} matrix (from \code{\link{build_DD}}).}
\item{d}{The matrix containing deterministic variable(s).}
\item{n_vars}{The number of variables.}
\item{n_lags}{The number of lags.}
\item{n_determ}{The number of deterministic terms.}
\item{n_reps}{The number of replications.}
}
\value{
\item{eval}{The evaluated density.}
}


