
\name{posteriors}
\alias{Z_posterior}
\alias{pi_sigma_posterior}
\alias{psi_posterior}
\title{Functions for drawing from the conditional posteriors}
\description{
The functions produce draws from the conditional posteriors and make up the blocks of the Gibbs sampler.
}
\usage{
pi_sigma_posterior(Z_r1, d, psi_r1, prior_pi, inv_prior_pi_omega, omega_pi, prior_s,
                   prior_nu, check_roots, n_vars, n_lags, n_T)
psi_posterior(Pi_r, Sigma_r, Z_r1, prior_psi, prior_psi_omega, D,
              n_vars, n_lags, n_determ)
Z_posterior(Y, d, Pi_r, Sigma_r, psi_r, Z_1, Lambda, n_vars, n_lags, n_T_, smooth_state)
}
\arguments{
\item{Pi_r}{The current draw of \code{Pi} (i.e. \code{Pi[,, r]}).}
\item{Sigma_r}{The current draw of \code{Sigma} (i.e. \code{Sigma[,, r]}).}
\item{psi_r1}{The previous draw of \code{psi} (i.e. \code{psi[r-1,]}).}
\item{psi_r}{The current draw of \code{psi} (i.e. \code{psi[r-1,]}).}
\item{Z_r1}{The previous draw of \code{Z} (i.e. \code{Z[,, r-1]}).}
\item{prior_pi}{The prior mean matrix for Pi.}
\item{inv_prior_pi_omega}{The inverse of the prior covariance matrix for Pi.}
\item{omega_pi}{The \code{inv_prior_pi_omega \%*\% prior_pi} matrix.}
\item{prior_s}{The prior for Sigma.}
\item{prior_nu}{The prior degrees of freedom.}
\item{prior_psi}{The prior mean for psi.}
\item{prior_psi_omega}{The prior covariance matrix for psi.}
\item{Y}{The data.}
\item{d}{The matrix containing deterministic variable(s).}
\item{D}{The \code{D} matrix (from \code{\link{build_DD}}).}
\item{Z_1}{The matrix \code{Z[1:n_lags,, 1]} (used as initial value).}
\item{Lambda}{The \code{Lambda} selection matrix.}
\item{check_roots}{If roots are to be checked to ensure stability.}
\item{smooth_state}{If the smoothed state should be computed and stored.}
\item{n_vars}{The number of variables.}
\item{n_lags}{The number of lags.}
\item{n_determ}{The number of deterministic terms.}
\item{n_T}{The number of time points in the original data (i.e. \code{nrow(Y)}).}
\item{n_T_}{The number of time points in the original data (i.e. \code{nrow(Y)}), minus the lag length.}
}

\value{
\code{pi_sigma_posterior} returns a list with:
\item{Pi_r}{The draw of \code{Pi}.}
\item{Sigma_r}{The draw of \code{Sigma}.}
\item{num_try}{The try at which a stable draw was obtained.}
\item{root}{The maximum eigenvalue (in modulus) of the system.}
\code{psi_posterior} returns:
\item{psi_r}{The draw of \code{psi}.}
\code{Z_posterior} returns a list with:
\item{Z_r}{The draw of \code{Z}.}
\item{smoothed_Z_r}{(Only if \code{smooth_state == TRUE}) The smoothed state.}
}
