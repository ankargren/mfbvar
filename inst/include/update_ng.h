#ifndef _UPDATE_NG_H_
#define _UPDATE_NG_H_

double posterior_phi_mu(const double lambda, const double phi_mu, const arma::vec omega,
                        const arma::uword nm) {
  double log_prob = arma::accu((phi_mu-1)*arma::log(omega)-0.5*lambda*phi_mu*omega) +
    nm*(phi_mu*std::log(lambda*phi_mu*0.5) -
    std::lgamma(phi_mu)) -
    exp(-phi_mu);

  return log_prob;
}

void update_ng(double & phi_mu, double & lambda_mu, arma::vec & omega, arma::uword nm,
               const double c0, const double c1, double s,
               const arma::vec & psi_i, const arma::vec & prior_psi_mean) {

  // Update omega
  double gig_lambda = phi_mu-0.5;
  double gig_chi = lambda_mu * phi_mu;
  arma::vec gig_psi = arma::pow(psi_i-prior_psi_mean, 2.0);
  for (arma::uword i = 0; i < nm; ++i) {
    omega(i) = rgig(gig_lambda, gig_chi, gig_psi(i));
  }

  // Update lambda
  lambda_mu = R::rgamma((double)nm * phi_mu + c0, 1.0/(0.5 * phi_mu * arma::accu(omega) + c1)); // Check parametrization

  // Update phi
  double phi_mu_proposal = phi_mu * std::exp(R::rnorm(0.0, s));
  double prob = exp(posterior_phi_mu(lambda_mu, phi_mu_proposal, omega, nm)-posterior_phi_mu(lambda_mu, phi_mu, omega, nm))*phi_mu_proposal/phi_mu;
  double u = R::runif(0.0, 1.0);
  if (u < prob) {
    phi_mu = phi_mu_proposal;
  }
}

#endif
