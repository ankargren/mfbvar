#include "mfbvar.h"

// [[Rcpp::export]]
double posterior_phi_mu(const double lambda, const double phi_mu, const arma::vec omega,
                        const arma::uword nm) {
  double log_prob = arma::accu((phi_mu-1)*arma::log(omega)-0.5*lambda*phi_mu*omega) +
    nm*(phi_mu*std::log(lambda*phi_mu*0.5) -
    std::lgamma(phi_mu)) - phi_mu;

  return log_prob;
}

void update_ng(double & phi_mu, double & lambda_mu, arma::vec & omega, arma::uword nm,
               const double c0, const double c1, double s,
               const arma::vec & psi_i, const arma::vec & prior_psi_mean, double & accept) {
  // phi_mu: the shrinkage parameter phi_mu
  // lambda_mu: the shrinkage parameter lambda_mu
  // omega: the idiosyncratic shrinkage parameters omega
  // nm: n_vars * n_determ (number of parameters)
  // c0: hyperparameter c0
  // c1: hyperparameter c1
  // s: scale of proposal
  // psi_i: the steady-state parameters
  // prior_psi_mean: the prior means of psi_i
  // accept: indicator for whether the proposal is accepted or not

  // Update omega
  double gig_lambda = phi_mu-0.5;
  //double gig_chi = lambda_mu * phi_mu;
  //arma::vec gig_psi = arma::pow(psi_i-prior_psi_mean, 2.0);

  arma::vec gig_chi = arma::pow(psi_i-prior_psi_mean, 2.0);
  double gig_psi = lambda_mu * phi_mu;

  for (arma::uword i = 0; i < nm; ++i) {
    //omega(i) = do_rgig1(gig_lambda, gig_chi, gig_psi(i));
    omega(i) = do_rgig1(gig_lambda, gig_chi(i), gig_psi);
  }
  // Update lambda
  lambda_mu = R::rgamma((double)nm * phi_mu + c0, 1.0/(0.5 * phi_mu * arma::accu(omega) + c1));

  // Update phi
  double phi_mu_proposal = phi_mu * std::exp(R::rnorm(0.0, s));
  double prob = exp(posterior_phi_mu(lambda_mu, phi_mu_proposal, omega, nm)-posterior_phi_mu(lambda_mu, phi_mu, omega, nm))*phi_mu_proposal/phi_mu;
  double u = R::runif(0.0, 1.0);
  if (u < prob) {
    phi_mu = phi_mu_proposal;
    accept = 1.0;
  } else {
    accept = 0.0;
  }
}

void update_s(double & s,
                arma::running_stat<double> & stats,
                double accept,
                arma::uword i,
                double M) {
  double s_prop;
  arma::vec min_vec(2);
  min_vec(0) = 0.01;
  double batch;
  double target = 0.44;
  stats(accept);
  if (i % 100 == 0) {
    batch = i/100;
    min_vec(1) = std::pow(batch, -0.5);
    if (stats.mean() > target) {
      s_prop = log(s) + arma::min(min_vec);
      if (s_prop < M){
        s = std::exp(s_prop);
      }
    } else {
      s_prop = log(s) - arma::min(min_vec);
      if (s_prop > -M){
        s = std::exp(s_prop);
      }
    }
    stats.reset();
  }
}
