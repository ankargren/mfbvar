#include "mfbvar.h"
#include "ss_utils.h"

arma::mat sample_Pi_mat(const arma::mat & post_Pi,
                        const arma::mat & post_Pi_Omega,
                        const arma::mat & Sigma,
                        bool check_roots,
                        arma::uword n_vars) {
  bool stationarity_check = false;
  int num_try = 0, iter = 0;
  double root = 1000;
  arma::mat Pi;
  arma::uword n_dim = post_Pi.n_rows;
  arma::mat Pi_comp = arma::mat(n_dim, n_dim, arma::fill::zeros);
  Pi_comp.submat(n_vars, 0, n_dim - 1, n_dim - n_vars - 1) =
    arma::eye(n_dim - n_vars, n_dim - n_vars);

  while (stationarity_check == false) {
    iter += 1;
    Pi = rmatn(post_Pi.t(), post_Pi_Omega, Sigma);
    if (check_roots) {
      Pi_comp.rows(0, n_vars-1) = Pi;
      root = max_eig_cpp(Pi_comp);
    } else {
      root = 0.0;
    }
    if (root < 1.0) {
      stationarity_check = true;
      num_try = iter;
    }
    if (iter == 1000) {
      Rcpp::stop("Attemped to draw stationary Pi 1,000 times.");
    }
  }
  return Pi;
}

arma::mat sample_Pi_vec(const arma::mat & Sigma_inv,
                        const arma::mat & X,
                        const arma::mat & y,
                        const arma::vec & prior_Pi_Omega_vec_inv,
                        const arma::mat & Omega_Pi,
                        bool check_roots,
                        arma::uword n_vars,
                        arma::uword n_lags) {
  arma::mat post_Pi_Omega_inv = arma::kron(Sigma_inv, X.t() * X);
  post_Pi_Omega_inv.diag() += prior_Pi_Omega_vec_inv;
  arma::mat L = arma::chol(post_Pi_Omega_inv, "lower");
  arma::vec b = arma::vectorise(X.t() * y * Sigma_inv + Omega_Pi);
  arma::vec u1 = arma::solve(arma::trimatl(L), b);
  arma::vec u2 = arma::solve(arma::trimatu(L.t()), u1);
  arma::mat u3 = arma::vec(n_vars*(n_vars*n_lags + 1));
  arma::mat Pi;
  arma::uword n_dim = n_vars * n_lags;
  arma::mat Pi_comp = arma::mat(n_dim, n_dim, arma::fill::zeros);
  Pi_comp.submat(n_vars, 0, n_dim - 1, n_dim - n_vars - 1) =
    arma::eye(n_dim - n_vars, n_dim - n_vars);
  bool stationarity_check = false;
  int num_try = 0, iter = 0;
  double root = 1000;
  while (stationarity_check == false) {
    iter += 1;
    u3.imbue(norm_rand);
    arma::vec u4 = arma::solve(arma::trimatu(L.t()), u3);
    arma::vec Pi_vec = u2 + u4;
    Pi = arma::trans(arma::reshape(Pi_vec, n_vars*n_lags, n_vars));
    if (check_roots) {
      Pi_comp.rows(0, n_vars-1) = Pi;
      root = max_eig_cpp(Pi_comp);
    } else {
      root = 0.0;
    }
    if (root < 1.0) {
      stationarity_check = true;
      num_try = iter;
    }
    if (iter == 1000) {
      Rcpp::stop("Attemped to draw stationary Pi 1,000 times.");
    }
  }
  return Pi;
}
