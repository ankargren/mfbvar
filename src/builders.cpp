// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "inst/include/mfbvarConfig.h"
using namespace Rcpp;
using namespace arma;

//' @describeIn build_U Build the U matrix (C++ implementation)
//' @templateVar n_vars TRUE
//' @templateVar n_lags TRUE
//' @keywords internal
//' @template man_template
// [[Rcpp::export]]
arma::mat build_U_cpp(arma::mat Pi, int n_determ, int n_vars, int n_lags){
  arma::mat U(n_vars*n_determ*(n_lags+1), n_vars*n_determ, fill::zeros);
  int pm = n_determ*n_vars;
  for(int i = 0; i < pm; i++){
    U(i,i) = 1;
  }
  arma::mat idmat(n_determ, n_determ, fill::eye);

  for(int i = 0; i < n_lags; i++){
    U.rows((i+1)*pm, (i+2)*pm - 1) = arma::kron(idmat, Pi.cols(i * n_vars, (i+1)*n_vars - 1));
  }

  return(U);

}


