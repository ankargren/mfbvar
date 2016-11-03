// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat build_U_cpp(mat Pi, int n_determ, int n_vars, int n_lags){
  mat U(n_vars*n_determ*(n_lags+1), n_vars*n_determ, fill::zeros);
  int pm = n_determ*n_vars;
  for(int i = 0; i < pm; i++){
    U(i,i) = 1;
  }
  mat idmat(n_determ, n_determ, fill::eye);

  for(int i = 0; i < n_lags; i++){
    U.rows((i+1)*pm, (i+2)*pm - 1) = kron(idmat, Pi.cols(i * n_vars, (i+1)*n_vars - 1));
  }

  return(U);

}


// [[Rcpp::export]]
mat build_demeaned_z_cpp(mat z, mat psi, mat d, int n_T, int n_vars){
  /* psi must be a row vector*/
  mat demeaned_z = z;
  mat idmat(n_T, n_T, fill::eye);
  mat idmat2(n_vars, n_vars, fill::eye);
  demeaned_z = z - kron(idmat, psi) * kron(vectorise(d.t()), idmat2);
  return(demeaned_z);
}

