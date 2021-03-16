#ifndef _UPDATE_IW_H_
#define _UPDATE_IW_H_
void update_iw(arma::mat & post_Pi_Omega,
               arma::mat & post_Pi,
               arma::mat & post_S,
               const arma::mat & X, const arma::mat & y,
               const arma::mat& prior_Pi_mean, const arma::mat & prior_Pi_Omega,
               const arma::mat & inv_prior_Pi_Omega,
               const arma::mat Omega_Pi, arma::mat prior_S);
#endif
