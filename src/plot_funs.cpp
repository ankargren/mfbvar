//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::export]]
void variances_fsv(arma::cube & variances, const arma::cube & latent, const arma::cube & facload, arma::uvec variables_num, arma::uword n_fac, arma::uword n_reps, arma::uword n_T, arma::uword n_vars, arma::uword n_plotvars) {
  arma::mat facload_i, fac_i, idi_i, variance_i;

  for (arma::uword i = 0; i < n_reps; ++i) {
    for (arma::uword tt = 0; tt < n_T; ++tt) {
      facload_i = facload.slice(i).rows(variables_num-1);
      fac_i = latent.slice(i).row(tt).cols(n_vars, n_vars+n_fac-1);
      idi_i = latent.slice(i).row(tt);
      idi_i = idi_i.cols(variables_num-1);
      variance_i = facload_i * arma::diagmat(arma::exp(fac_i)) * facload_i.t();
      variance_i.diag() += arma::exp(idi_i);
      variances.slice(i).row(tt) = arma::sqrt(variance_i.diag().t());
    }
  }
}

//[[Rcpp::export]]
void variances_csv(arma::cube & variances, const arma::cube & Sigma, const arma::mat & f, arma::uword n_T, arma::uword n_reps, arma::uvec variables_num) {
  arma::vec Sigma_i;
  for (arma::uword i = 0; i < n_T; ++i) {
    for (arma::uword j = 0; j < n_reps; ++j) {
      Sigma_i = Sigma.slice(j).diag();
      variances.slice(j).row(i) = std::exp(0.5 * f(j, i)) * arma::sqrt(Sigma_i.rows(variables_num-1).t());
    }
  }
}
