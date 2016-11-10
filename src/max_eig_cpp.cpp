#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double max_eig_cpp(arma::mat M) {
  arma::cx_vec eigval = arma::eig_gen(M);
  return abs(max(eigval));
}
