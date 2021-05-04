#include "mfbvar.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @title Find maximum eigenvalue
//'
//' @description The function computes the maximum eigenvalue.
//' @aliases max_eig max_eig_cpp
//' @param A Symmetrix matrix whose maximum eigenvalue is to be computed
//' @keywords internal
//' @noRd
//' @return The maximum eigenvalue.
// [[Rcpp::export]]
double max_eig_cpp(const arma::mat & A) {
  arma::cx_vec eigval = arma::eig_gen(A);
  return max(abs(eigval));
}
