\name{prior_pi_sigma}
\alias{prior_pi_sigma}
\title{Create the priors for Pi and Sigma}
\description{
Creates the prior mean and covariance for Pi given the hyperparameters, and the prior parameters for Sigma.
}
\usage{
prior_pi_sigma(lambda1, lambda2, prior_mean, Y, n_lags, nu)
}
\arguments{
\item{lambda1}{The overall tightness.}
\item{lambda2}{The lag decay.}
\item{prior_mean}{The prior means for the AR(1) coefficients.}
\item{Y}{The data.}
\item{n_lags}{The number of lags.}
\item{nu}{The parameter \eqn{\nu}.}
}
\value{
\item{prior_pi}{The prior mean matrix for Pi.}
\item{prior_pi_omega}{The prior covariance matrix for Pi.}
\item{prior_s}{The prior for Sigma.}
}

