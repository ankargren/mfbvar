#include "mfbvar.h"

//' @rdname dnorminvwish
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rmatn(const arma::mat & M, const arma::mat & Q, const arma::mat & P){
/*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int p = P.n_rows;
  int q = Q.n_rows;
  arma::mat L = arma::chol(Q, "upper");
  arma::mat C = arma::chol(P, "lower");
  arma::mat X = arma::reshape(arma::vec(Rcpp::rnorm(p * q)), p, q);
  X = M + C * X * L;
  return(X);
}


//' @rdname dnorminvwish
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rinvwish(int v, const arma::mat & S){
  int p = S.n_rows;
  arma::mat L = arma::chol(arma::inv_sympd(S), "lower");
  arma::mat A(p,p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df));
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  arma::mat LA_inv = arma::inv(arma::trimatl(arma::trimatl(L) * arma::trimatl(A)));
  arma::mat X = LA_inv.t() * LA_inv;

  return(X);
}


//' @rdname dmultn
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rmultn(const arma::vec & m, const arma::mat & Sigma){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int p = Sigma.n_rows;
  arma::vec X = Rcpp::rnorm(p);
  arma::mat L = arma::chol(Sigma, "lower");
  X = m + L * X;
  return(X);
}



