#include "mfbvar.h"
void update_dl(arma::mat & prior_Pi_Omega, arma::vec & aux,
               arma::vec & local, double & global, const arma::mat & Pi_i,
               arma::uword n_vars, arma::uword n_lags, const double a,
               arma::vec & slice, bool fixate_aux, bool fixate_global,
               bool fixate_local, bool gig = true, bool intercept = true) {

  arma::vec Pi_vec;
  if (intercept) {
    Pi_vec = arma::vectorise(Pi_i.rows(1, n_vars*n_lags));
  } else {
    Pi_vec = arma::vectorise(Pi_i);
  }


  arma::uword K = Pi_vec.n_elem;
  if (!fixate_aux) {
    for (arma::uword i = 0; i < K; ++i) {
      aux[i] = 1.0/rig(global * local[i] / fabs(Pi_vec[i]), 1.0);
    }
  }
  arma::vec Pi_local = arma::abs(Pi_vec) / local;

  if (!fixate_global) {
    global = do_rgig1(K*(a-1.0), 2.0 * arma::accu(Pi_local), 1.0);
  }

  if (!fixate_local) {
    if (gig) {
      for (arma::uword i = 0; i < K; ++i) {
        local[i] = do_rgig1((a-1.0), 2.0 * fabs(Pi_vec[i]), 1.0);
      }
    } else {
      arma::vec u1 = arma::vec(K);
      std::generate(u1.begin(), u1.end(), ::unif_rand);
      u1 %= arma::exp(-0.5 / slice);
      arma::vec lb = 0.5/(arma::log(1/u1));
      double Flb;
      arma::vec u2 = arma::vec(K);
      for (arma::uword i = 0; i < K; ++i) {
        Flb = R::pgamma(lb[i], 1-a, 1/fabs(Pi_vec[i]), true, false);
        u2[i] = R::runif(Flb, 1.0);
      }
      arma::uvec u3 = arma::find(u2 > 1-(1e-16));
      if (u3.n_elem > 0) {
        u2.elem(u3).fill(1-(1e-16));
      }
      for (arma::uword i = 0; i < K; ++i) {
        slice[i] = R::qgamma(u2[i], 1-a, 1/fabs(Pi_vec[i]), true, false);
      }
      local = 1/slice;
    }

    local = local / arma::accu(local);

    arma::uvec local_idx = arma::find(local < 1e-20);
    local.elem(local_idx).fill(1e-20);
  }



  if (intercept) {
    prior_Pi_Omega.rows(1, n_vars*n_lags) = arma::reshape(aux % arma::pow(global * local, 2.0), n_vars*n_lags, n_vars);
  } else {
    prior_Pi_Omega = arma::reshape(aux % arma::pow(global * local, 2.0), n_vars*n_lags, n_vars);
  }


}
