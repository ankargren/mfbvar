// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat rmatnorm(mat M, mat Q, mat P){
/*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  RNGScope scope;
  int p = P.n_rows;
  int q = Q.n_rows;
  mat L = chol(Q, "upper");
  mat C = chol(P, "lower");
  mat X = reshape(vec(rnorm(p * q)), p, q);
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
mat rinvwish(int v, mat S){
  RNGScope scope;
  int p = S.n_rows;
  mat L = chol(inv_sympd(S), "lower");
  mat A(p,p, fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df));
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
  mat X = LA_inv.t() * LA_inv;

  return(X);
}


// [[Rcpp::export]]
vec rmultn(vec m, mat Sigma){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
   RNGScope scope;
  int p = Sigma.n_rows;
  vec X = rnorm(p);
  mat L = chol(Sigma, "lower");
  X = m + L * X;
  return(X);
}
