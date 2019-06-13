#include "mfbvar.h"
void update_dl(arma::mat & prior_Pi_Omega, arma::vec & aux,
               arma::vec & local, double & global, const arma::mat & Pi_i,
               arma::uword n_vars, arma::uword n_lags, const double a,
               bool intercept = true) {

  arma::vec Pi_vec;
  if (intercept) {
    Pi_vec = arma::vectorise(Pi_i.rows(1, n_vars*n_lags));
  } else {
    Pi_vec = arma::vectorise(Pi_i);
  }


  arma::uword K = Pi_vec.n_elem;

  for (arma::uword i = 0; i < K; ++i) {
    aux[i] = 1.0/rig(global * local[i] / fabs(Pi_vec[i]), 1.0);
  }
  arma::vec Pi_local = arma::abs(Pi_vec) / local;

  global = do_rgig1(K*(a-1.0), 2.0 * arma::accu(Pi_local), 1.0);

  for (arma::uword i = 0; i < K; ++i) {
    local[i] = do_rgig1((a-1.0), 2.0 * fabs(Pi_vec[i]), 1.0);
  }

  local = local / arma::accu(local);
  if (intercept) {
    prior_Pi_Omega.rows(1, n_vars*n_lags) = arma::reshape(aux % arma::pow(global * local, 2.0), n_vars*n_lags, n_vars);
  } else {
    prior_Pi_Omega = arma::reshape(aux % arma::pow(global * local, 2.0), n_vars*n_lags, n_vars);
  }


}
