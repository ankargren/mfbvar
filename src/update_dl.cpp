#include "mfbvar.h"
void update_dl(arma::mat & prior_Pi_Omega, arma::vec & aux,
               arma::vec & local, double & global, const arma::mat & Pi_i,
               arma::uword n_vars, arma::uword n_lags, const double a) {


  arma::vec Pi_vec = arma::vectorise(Pi_i.rows(1, n_vars*n_lags));
  arma::uword K = Pi_vec.n_elem;
  arma::vec local_Pi = local / arma::abs(Pi_vec);

  for (arma::uword i = 0; i < K; ++i) {
    aux[i] = 1.0/rig(global * local_Pi[i], 1.0);
  }

  global = do_rgig1(K*(a-1.0), 1.0, 2.0 * arma::accu(1.0/local_Pi));

  for (arma::uword i = 0; i < K; ++i) {
    local[i] = do_rgig1((a-1.0), 1.0, 2.0 * fabs(Pi_vec[i]));
  }
  local = local * (1.0 / arma::accu(local));
  prior_Pi_Omega.rows(1, n_vars*n_lags) = arma::reshape(aux * arma::pow(global * local, 2.0), n_vars*n_lags, n_vars);

}
