#include "mfbvar.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @title Find maximum eigenvalue
//'
//' @description The function computes the maximum eigenvalue.
//' @aliases max_eig max_eig_cpp
//' @templateVar A TRUE
//' @template man_template
//' @keywords internal
//' @noRd
//' @return The maximum eigenvalue.
// [[Rcpp::export]]
double max_eig_cpp(const arma::mat & A) {
  arma::cx_vec eigval = arma::eig_gen(A);
  return max(abs(eigval));
}
