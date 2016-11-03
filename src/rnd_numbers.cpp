// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat rmatn(arma::mat M, arma::mat Q, arma::mat P){
/*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  RNGScope scope;
  int p = P.n_rows;
  int q = Q.n_rows;
  arma::mat L = chol(Q, "upper");
  arma::mat C = chol(P, "lower");
  arma::mat X = reshape(arma::vec(rnorm(p * q)), p, q);
  X = M + C * X * L;
  return(X);
}


/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition
#--------------------------------------------------------
# NOTE: output is identical to riwish from MCMCpack
#       provided the same random seed is used
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix
#   v     degrees of freedom
#-------------------------------------------------------*/
// [[Rcpp::export]]
arma::mat rinvwish(int v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = arma::chol(arma::inv_sympd(S), "lower");
  arma::mat A(p,p, fill::zeros);
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


// [[Rcpp::export]]
arma::vec rmultn(arma::vec m, arma::mat Sigma){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
   RNGScope scope;
  int p = Sigma.n_rows;
  arma::vec X = rnorm(p);
  arma::mat L = arma::chol(Sigma, "lower");
  X = m + L * X;
  return(X);
}
