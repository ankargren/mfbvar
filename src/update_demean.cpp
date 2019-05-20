#include <RcppArmadillo.h>
// [[Rcpp::export]]
void update_demean(arma::mat & my, arma::mat & mu_long,
                   const arma::mat & y_in_p, const arma::mat & mu_mat, const arma::mat & d1,
                   const arma::mat & Psi_i, const arma::mat & Lambda_single,
                   arma::uword n_vars, arma::uword n_q, arma::uword n_Lambda, arma::uword n_T) {
  my.cols(0, n_vars - n_q - 1) = y_in_p.cols(0, n_vars - n_q - 1) - mu_mat.cols(0, n_vars - n_q - 1);
  mu_long.rows(0, n_Lambda-1) = d1.tail_rows(n_Lambda) * Psi_i.t();
  mu_long.rows(n_Lambda, n_T+n_Lambda-1) = mu_mat;
  for (arma::uword j = 0; j < n_T; ++j) {
    my.row(j).cols(n_vars - n_q - 1, n_vars - 1) = y_in_p.row(j).cols(n_vars - n_q - 1, n_vars - 1) - Lambda_single * mu_long.rows(j, j+n_Lambda-1).cols(n_vars - n_q - 1, n_vars - 1);// Needs fixing
  }
}
