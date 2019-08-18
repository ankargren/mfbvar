#include "mfbvar.h"

// [[Rcpp::export]]
arma::mat rsimsm_adaptive_univariate(arma::mat y_, arma::mat Phi, arma::mat Sigma,
                                         arma::mat Lambda, arma::mat Z1, arma::uword n_q_, arma::uword T_b, arma::mat f) {
  return simsm_adaptive_univariate(y_, Phi, Sigma, Lambda, Z1, n_q_, T_b, f);
}
