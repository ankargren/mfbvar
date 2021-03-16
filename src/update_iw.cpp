#include "mfbvar.h"
void update_iw(arma::mat & post_Pi_Omega,
               arma::mat & post_Pi,
               arma::mat & post_S,
               const arma::mat & X, const arma::mat & y,
               const arma::mat& prior_Pi_mean, const arma::mat & prior_Pi_Omega,
               const arma::mat & inv_prior_Pi_Omega,
               const arma::mat Omega_Pi, arma::mat prior_S) {
  arma::mat XX = X.t() * X;
  arma::mat XX_inv = arma::inv_sympd(XX);
  arma::mat Pi_sample = XX_inv * (X.t() * y);
  post_Pi_Omega = arma::inv_sympd(inv_prior_Pi_Omega + XX);
  post_Pi = post_Pi_Omega * (Omega_Pi + X.t() * y);
  arma::mat S = arma::trans((y - X * Pi_sample)) * (y - X * Pi_sample);
  arma::mat Pi_diff = prior_Pi_mean - Pi_sample;
  post_S = prior_S + S +
    Pi_diff.t() * arma::inv_sympd(prior_Pi_Omega + XX_inv) * Pi_diff;
}
