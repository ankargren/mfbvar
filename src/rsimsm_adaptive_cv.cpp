#include "mfbvar.h"

// [[Rcpp::export]]
arma::mat rsimsm_adaptive_cv(arma::mat y_, arma::mat Phi, arma::mat Sigma_chol,
                                         arma::mat Lambda, arma::mat Z1, arma::uword n_q_, arma::uword T_b) {
  return simsm_adaptive_cv(y_, Phi, Sigma_chol, Lambda, Z1, n_q_, T_b);
}

// [[Rcpp::export]]
arma::mat rsimsm_adaptive_sv(arma::mat y_, arma::mat Phi, arma::cube Sigma_chol,
                             arma::mat Lambda, arma::mat Z1, arma::uword n_q_, arma::uword T_b) {
  return simsm_adaptive_sv(y_, Phi, Sigma_chol, Lambda, Z1, n_q_, T_b);
}

