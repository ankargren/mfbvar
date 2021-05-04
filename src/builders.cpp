#include "mfbvar.h"

//' @describeIn build_U Build the U matrix (C++ implementation)
//' @param n_vars Number of variables
//' @param n_lags Number of lags
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags){
  arma::mat U(n_vars*n_determ*(n_lags+1), n_vars*n_determ, arma::fill::zeros);
  int pm = n_determ*n_vars;
  for(int i = 0; i < pm; i++){
    U(i,i) = 1;
  }

  for(int i = 0; i < n_lags; i++){
    for (int j = 0; j < n_determ; j++) {
      U(arma::span((i+1)*pm + j*n_vars,(i+1)*pm + (j+1)*n_vars-1), arma::span(j*n_vars, (j+1)*n_vars-1)) = Pi.cols(i * n_vars, (i+1)*n_vars - 1);
    }
  }

  return(U);

}

// [[Rcpp::export]]
arma::mat create_X(const arma::mat & y, arma::uword k) {
  arma::uword TT = y.n_rows - k;
  arma::uword n = y.n_cols;
  arma::mat X = arma::mat(TT, n*k + 1, arma::fill::ones);

  for (arma::uword i = 0; i < k; i++) {
    X.cols(i*n+1,(i+1)*n) = y.rows(k-i-1, TT + k - 2 - i);
  }
  return X;
}

// [[Rcpp::export]]
arma::mat create_X_noint(const arma::mat & y, arma::uword k) {
  arma::uword TT = y.n_rows - k;
  arma::uword n = y.n_cols;
  arma::mat X = arma::mat(TT, n*k, arma::fill::zeros);

  for (arma::uword i = 0; i < k; i++) {
    X.cols(i*n,(i+1)*n-1) = y.rows(k-i-1, TT + k - 2 - i);
  }
  return X;
}

// [[Rcpp::export]]
arma::mat create_X_t(const arma::mat & y) {
  arma::uword np = y.n_elem;
  arma::mat X = arma::mat(np+1, 1, arma::fill::ones);
  X.rows(1, np) = arma::reshape(arma::trans(arma::flipud(y)), np, 1);
  return X;
}

// [[Rcpp::export]]
arma::mat create_X_t_noint(const arma::mat & y) {
  arma::uword np = y.n_elem;
  arma::mat X = arma::reshape(arma::trans(arma::flipud(y)), np, 1);
  return X;
}
