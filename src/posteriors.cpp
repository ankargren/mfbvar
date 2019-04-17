#include "mfbvar.h"

// [[Rcpp::export]]
arma::mat posterior_psi_Omega_fsv(const arma::mat & U, const arma::mat & D_mat,
                                  const arma::mat & idivar, const arma::mat & inv_prior_psi_Omega) {
  arma::mat mid_mat = arma::mat(D_mat.n_cols * idivar.n_cols, D_mat.n_cols * idivar.n_cols, arma::fill::zeros);
  arma::uword n_T = D_mat.n_rows;
  arma::mat D_temp;
  for (arma::uword i = 0; i < n_T; i++) {
    D_temp = D_mat.row(i);
    mid_mat += arma::kron(D_temp.t() * D_temp, arma::diagmat(arma::pow(idivar.row(i), -1.0)));
  }
  arma::mat psi_Omega = arma::inv_sympd((U.t() * mid_mat) * U + inv_prior_psi_Omega);
  return psi_Omega;
}



// [[Rcpp::export]]
arma::vec posterior_psi_mean_fsv(const arma::mat & U, const arma::mat & D_mat, const arma::mat & idivar,
                                 const arma::vec & inv_prior_psi_Omega_mean, const arma::mat & post_psi_Omega,
                                 const arma::mat & Y_tilde) {
  arma::mat SigmaYD = arma::mat(idivar.n_cols, D_mat.n_cols, arma::fill::zeros);
  arma::uword n_T = D_mat.n_rows;
  for (arma::uword i = 0; i < n_T; i++) {
    SigmaYD += arma::trans(Y_tilde.row(i) / idivar.row(i)) * D_mat.row(i);
  }
  arma::vec sigmaYD = arma::vectorise(SigmaYD);
  arma::vec psi = post_psi_Omega * (U.t() * sigmaYD + inv_prior_psi_Omega_mean);
  return psi;
}


// [[Rcpp::export]]
void posterior_psi_fsv(arma::vec & psi_i, arma::mat & mu_mat,
                             const arma::mat & Pi_i, const arma::mat & D_mat,
                             const arma::mat & idivar, const arma::mat & inv_prior_psi_Omega,
                             const arma::mat & Z_i, const arma::mat & X,
                             const arma::mat & startfacload, const arma::mat & startfac,
                             const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                             int n_determ, int n_vars, int n_lags) {
  arma::mat U = build_U_cpp(Pi_i, n_determ, n_vars, n_lags);
  arma::mat post_psi_Omega = posterior_psi_Omega_fsv(U, D_mat, idivar, inv_prior_psi_Omega);
  arma::mat Y_tilde = Z_i - X * Pi_i.t() - arma::trans(startfacload * startfac);

  arma::mat post_psi = posterior_psi_mean_fsv(U, D_mat, idivar, inv_prior_psi_Omega_mean,
                                              post_psi_Omega, Y_tilde);
  psi_i = rmultn(post_psi, post_psi_Omega);
  arma::mat Psi_i = arma::mat(psi_i.begin(), n_vars, n_determ);

  mu_mat = dt * Psi_i.t();
}


// [[Rcpp::export]]
arma::vec posterior_psi_mean_iw(const arma::mat & U, const arma::mat & D_mat, const arma::mat & Sigma_i,
                                const arma::vec & inv_prior_psi_Omega_mean, const arma::mat & post_psi_Omega,
                                const arma::mat & Y_tilde) {

  arma::mat SigmaYD = arma::inv_sympd(Sigma_i) * (Y_tilde.t() * D_mat);
  arma::vec sigmaYD = arma::vectorise(SigmaYD);
  arma::vec psi = post_psi_Omega * (U.t() * sigmaYD + inv_prior_psi_Omega_mean);
  return psi;
}

// [[Rcpp::export]]
arma::mat posterior_psi_Omega_iw(const arma::mat & U, const arma::mat & D_mat,
                                 const arma::mat & Sigma_i, const arma::mat & inv_prior_psi_Omega) {
  arma::mat mid_mat = arma::kron(D_mat.t() * D_mat, arma::inv_sympd(Sigma_i));
  arma::mat psi_Omega = arma::inv_sympd((U.t() * mid_mat) * U + inv_prior_psi_Omega);
  return psi_Omega;
}

// [[Rcpp::export]]
void posterior_psi_iw(arma::vec & psi_i, arma::mat & mu_mat,
                           const arma::mat & Pi_i, const arma::mat & D_mat,
                           const arma::mat & Sigma_i, const arma::mat & inv_prior_psi_Omega,
                           const arma::mat & Z_i, const arma::mat & X,
                           const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                           int n_determ, int n_vars, int n_lags) {
  arma::mat U = build_U_cpp(Pi_i, n_determ, n_vars, n_lags);
  arma::mat post_psi_Omega = posterior_psi_Omega_iw(U, D_mat, Sigma_i, inv_prior_psi_Omega);
  arma::mat Y_tilde = Z_i - X * Pi_i.t();
  arma::mat post_psi = posterior_psi_mean_iw(U, D_mat, Sigma_i, inv_prior_psi_Omega_mean,
                                             post_psi_Omega, Y_tilde);
  psi_i = rmultn(post_psi, post_psi_Omega);
  arma::mat Psi_i = arma::mat(psi_i.begin(), n_vars, n_determ);

  mu_mat = dt * Psi_i.t();
}
