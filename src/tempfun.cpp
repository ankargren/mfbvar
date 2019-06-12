#include "mfbvar.h"
#include "minn_utils.h"
#include "update_dl.h"
#include "mvn_par.h"
// [[Rcpp::export]]
void mcmc_minn_fsv_tmp(const arma::mat & y_in,
                   arma::cube& Pi, arma::mat armah,
                   arma::mat & aux, arma::vec & global, arma::mat & local,
                   arma::mat prior_Pi_Omega,
                   const arma::vec& prior_Pi_AR1,
                   arma::uword n_fac, arma::uword n_reps,
                   arma::uword n_lags, arma::uword n_vars,
                   arma::uword n_T,
                   const double a) {

  arma::mat Pi_i = Pi.slice(0);
  arma::mat X;

  arma::mat y_i = y_in.rows(n_lags, y_in.n_rows-1);

  arma::mat eps = arma::mat(n_vars*n_lags+1, n_vars);
  X = create_X(y_in, n_lags);

  double global_i = global(0);
  arma::vec aux_i = aux.row(0).t();
  arma::vec local_i = local.row(0).t();

  for (arma::uword i = 0; i < n_reps; ++i) {
    eps.imbue(norm_rand);
    arma::mat output(n_vars*n_lags+1, n_vars);
    Pi_parallel_rue Pi_parallel_i(output, y_i, X, prior_Pi_Omega, eps,
                                  armah, prior_Pi_AR1, n_T, n_vars, n_lags);
    RcppParallel::parallelFor(0, n_vars, Pi_parallel_i);

    Pi_i = output.t();
    Pi.slice(i) = Pi_i;
    update_dl(prior_Pi_Omega, aux_i, local_i, global_i, Pi_i.t(), n_vars, n_lags, a);
    global(i) = global_i;
    aux.row(i) = aux_i.t();
    local.row(i) = local_i.t();

  }

}
