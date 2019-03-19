#include "mfbvar.h"

//' @describeIn build_U Build the U matrix (C++ implementation)
//' @templateVar n_vars TRUE
//' @templateVar n_lags TRUE
//' @keywords internal
//' @template man_template
// [[Rcpp::export]]
arma::mat build_U_cpp(arma::mat Pi, int n_determ, int n_vars, int n_lags){
  arma::mat U(n_vars*n_determ*(n_lags+1), n_vars*n_determ, arma::fill::zeros);
  int pm = n_determ*n_vars;
  for(int i = 0; i < pm; i++){
    U(i,i) = 1;
  }
  arma::mat idmat(n_determ, n_determ, arma::fill::eye);

  for(int i = 0; i < n_lags; i++){
    U.rows((i+1)*pm, (i+2)*pm - 1) = arma::kron(idmat, Pi.cols(i * n_vars, (i+1)*n_vars - 1));
  }

  return(U);

}

// [[Rcpp::export]]
arma::mat create_X(arma::mat y, arma::uword k) {
  arma::uword TT = y.n_rows - k;
  arma::uword n = y.n_cols;
  arma::mat X = arma::mat(TT, n*k + 1, arma::fill::ones);

  for (arma::uword i = 0; i < k; i++) {
    X.cols(i*n+1,(i+1)*n) = y.rows(k-i-1, TT + k - 2 - i);
  }
  return X;
}

