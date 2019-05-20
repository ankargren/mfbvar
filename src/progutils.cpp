// Copyright of original code (excl. 'rtruncnorm'): Gregor Kastner (stochvol package)
// Copyright of modified code: Sebastian Ankargren (mfbvar package)
// The following code is a derivative work of the code
// developed by Gregor Kastner for the stochvol package, which
// is licensed GPL>=2. This code is therefore licensed under
// the terms of the GNU Public License, version 3.
#include <RcppArmadillo.h>
#include "progutils.h"
// b)
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(
    const arma::vec& omega_diag,
    double omega_offdiag,
    arma::vec& chol_diag,
    arma::vec& chol_offdiag) {
 chol_diag[0] = std::pow(omega_diag[0], 0.5);  // maybe speed up via iterators?
 for (int j = 1; j < int(omega_diag.size()); j++) {
  chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
  chol_diag[j] = std::pow(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1], 0.5);
 }
}

// Solves Chol*x = covector ("forward algorithm")
// [[Rcpp::export]]
void forwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector,
    arma::vec& htmp) {
 htmp[0] = covector[0]/chol_diag[0];
 for (int j = 1; j < int(chol_diag.size()); j++) {
  htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
 }
}

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp,
    arma::vec& h) {
 int T = chol_diag.size();
 h[T-1] = htmp[T-1]/chol_diag[T-1];
 for (int j = T-2; j >= 0; j--) {
  h[j] = (htmp[j] - chol_offdiag[j]*h[j+1])/chol_diag[j];
 }
}

double rtruncnorm(double m, double v) {
  double proposal = 1000;
  int iter = 0;
  while(std::abs(proposal)>1.0) {
    proposal = R::rnorm(m, std::pow(v, 0.5));
    iter += 1;
    if (iter > 1000) {
      Rcpp::stop("Unable to draw stationary phi.");
    }
  }
  return proposal;
}
