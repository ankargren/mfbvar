#include "mfbvar.h"
#include "minn_utils.h"
void fcst_csv(arma::mat & Z_fcst, const arma::mat & Z, const arma::mat & Pi,
              const arma::mat & Sigma_chol, double phi, double sigma, double vol_pred,
              arma::uword n_fcst, arma::uword n_lags, arma::uword n_vars) {
  arma::vec errors = arma::vec(n_vars);
  arma::mat x;
  Z_fcst.head_cols(n_lags) = Z.tail_rows(n_lags).t();
  for (arma::uword h = 0; h < n_fcst; ++h) {
    vol_pred = phi * vol_pred + R::rnorm(0.0, sigma);
    errors.imbue(norm_rand);
    errors = errors * std::exp(0.5 * vol_pred);
    x = create_X_t(Z_fcst.cols(0+h, n_lags-1+h).t(), false);
    Z_fcst.col(n_lags + h) = Pi * x + Sigma_chol * errors;
  }
}
