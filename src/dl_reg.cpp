#include "mfbvar.h"
#include "update_dl.h"
// [[Rcpp::export]]
void dl_reg(const arma::mat & y, arma::mat & x, arma::mat & beta,
            arma::mat & aux, arma::vec & global, arma::mat & local,
            arma::mat & prior_Pi_Omega, arma::uword n_reps,
            const double a, bool gig) {

  arma::mat eps = arma::mat(x.n_cols, 1);
  arma::mat beta_i = beta.row(0).t();

  double global_i = global(0);
  arma::vec aux_i = aux.row(0).t();
  arma::vec local_i = local.row(0).t();
  arma::vec slice = arma::vec(local_i.n_elem).fill(1.0);

  arma::mat Sigma, Sigma_inv, L;
  arma::vec mu;
  arma::uword n_lags = x.n_cols;

  for (arma::uword i = 0; i < n_reps; ++i) {

    eps.imbue(norm_rand);
    Sigma = x.t() * x;
    Sigma.diag() += 1/prior_Pi_Omega;
    Sigma_inv = arma::inv_sympd(Sigma);
    L = arma::chol(Sigma_inv, "lower");
    mu = Sigma_inv * x.t() * y;
    beta_i.col(0) = mu + L * eps;

    //beta_i.col(0) = mvn_rue(x, prior_Pi_Omega, y);
    beta.row(i) = beta_i.col(0).t();
    update_dl(prior_Pi_Omega, aux_i, local_i, global_i, beta_i, 1, n_lags, a, slice,
              false, false, false, gig, false);
    global(i) = global_i;
    aux.row(i) = aux_i.t();
    local.row(i) = local_i.t();
  }

}
