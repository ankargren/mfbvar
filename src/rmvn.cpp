#include "mfbvar.h"
#include "mvn.h"

// [[Rcpp::export]]
arma::vec rmvn(const arma::mat & Phi, const arma::vec & d, const arma::vec & alpha) {
    arma::uword n = Phi.n_rows;
    arma::uword p = Phi.n_cols;
    arma::vec theta(p, arma::fill::zeros);

    if (p > 1.1*n) {
      theta = mvn_bcm(Phi, d, alpha);
    } else {
      theta = mvn_rue(Phi, d, alpha);
    }

    return theta;
}

// [[Rcpp::export]]
arma::vec rmvn_ccm(const arma::mat & Phi, const arma::vec & d, const arma::vec & alpha, double c, double j) {
  return mvn_ccm(Phi, d, alpha, c, j);
}
