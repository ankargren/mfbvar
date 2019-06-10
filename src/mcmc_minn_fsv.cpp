#include "mfbvar.h"
#include "minn_utils.h"
#include "update_fsv.h"
#include "mvn.h"
// [[Rcpp::export]]
void mcmc_minn_fsv(const arma::mat & y_in_p,
                   arma::cube& Pi, arma::cube& Z, arma::cube& Z_fcst,
                   arma::mat& mu, arma::mat& phi, arma::mat& sigma,
                   arma::cube& f, arma::cube& facload, arma::cube& h,
                   const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                   const arma::vec& prior_Pi_AR1, const arma::mat& Z_1,
                   double bmu, double Bmu, double a0idi, double b0idi, double a0fac, double b0fac,
                   const Rcpp::NumericVector & Bsigma, double B011inv, double B022inv,
                   const Rcpp::NumericVector & priorh0, const arma::imat & armarestr,
                   const arma::mat & armatau2, // armatau2 is the matrix with prior variance of factor loadings
                   arma::uword n_fac, arma::uword n_reps,
                   arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                   arma::uword n_T, arma::uword n_fcst, arma::uword n_thin, bool verbose) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps, verbose);

  arma::mat Pi_i = Pi.slice(0);
  arma::mat X;
  arma::mat y_i = y_in_p;
  arma::mat x;
  arma::vec vol_pred;


  // fsv
  Rcpp::NumericMatrix curpara = Rcpp::NumericMatrix(3, n_vars + n_fac);
  arma::mat curpara_arma(curpara.begin(), curpara.nrow(), curpara.ncol(), false);
  curpara_arma.fill(0.0);
  curpara_arma.row(0).cols(0, n_vars - 1) = mu.col(0).t();
  curpara_arma.row(1) = phi.col(0).t();
  curpara_arma.row(2) = sigma.col(0).t();

  arma::vec mu_i = mu.col(0);
  arma::vec phi_i = phi.col(0);
  arma::vec sigma_i = sigma.col(0);

  arma::mat armaf = f.slice(0);
  arma::mat armafacload = facload.slice(0);
  arma::mat armah = h.slice(0);
  arma::mat cc_i = armaf.t() * armafacload.t();

  arma::vec armah0 = arma::vec(n_vars + n_fac);

  arma::mat Sig_i, y_hat, latent_nofac, h_j, X_j, y_j;
  arma::vec error_pred;
  arma::vec errors_sv = arma::vec(n_vars + n_fac);
  arma::vec errors_var = arma::vec(n_vars + n_fac);

  arma::mat Z_i = arma::mat(n_lags + y_in_p.n_rows, n_vars, arma::fill::zeros);
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  Z_i.rows(0, n_lags - 1) = Z_1;

  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags);
  }

  for (arma::uword i = 0; i < n_reps; ++i) {
    if (!single_freq) {
      Sig_i = arma::exp(0.5 * armah.head_cols(n_vars));
      y_i = simsm_adaptive_univariate(y_in_p, Pi_i, Sig_i, Lambda_comp, Z_1, n_q, T_b, cc_i);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags);
    }

    y_hat = y_i - X * Pi_i.t();

    if (i % n_thin == 0) {

      mu_i = curpara_arma.row(0).t();
      phi_i = curpara_arma.row(1).t();
      sigma_i = curpara_arma.row(2).t(); // sigma, not sigma2
      if (n_fcst > 0) {
        vol_pred = armah.tail_rows(1).t();
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors_sv.imbue(norm_rand);
          errors_var.imbue(norm_rand);
          vol_pred = mu_i + phi_i % (vol_pred - mu_i) + sigma_i % errors_sv; // Twice because we first need it for the volatility, then for the VAR
          error_pred = arma::exp(0.5 * vol_pred) % errors_var;
          x = create_X_t(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + armafacload * error_pred.tail_rows(n_fac) + error_pred.head_rows(n_vars);
        }
        Z_fcst.slice(i/n_thin) = Z_fcst_i.t();
      }


      Z.slice(i/n_thin) = Z_i;

      Pi.slice(i/n_thin) = Pi_i;


      f.slice(i/n_thin) = armaf;
      facload.slice(i/n_thin) = armafacload;
      h.slice(i/n_thin) = armah;


      mu.col(i/n_thin) = mu_i.head(n_vars);
      phi.col(i/n_thin) = phi_i;
      sigma.col(i/n_thin) = sigma_i;
    }
    update_fsv(armafacload, armaf, armah, armah0, curpara, armatau2, y_hat.t(),
               bmu, Bmu, a0idi, b0idi, a0fac, b0fac, Bsigma, B011inv, B022inv,
               priorh0, armarestr);

    cc_i = armaf.t() * armafacload.t(); // Common component
    latent_nofac = y_i - cc_i;

    for (arma::uword j = 0; j < n_vars; j++) {
      arma::vec h_j = arma::exp(-0.5 * armah.col(j));
      arma::mat X_j = X.each_col() % h_j;
      arma::vec y_j = latent_nofac.col(j) % h_j;
      Pi_i.row(j) = arma::trans(mvn_ccm(X_j, prior_Pi_Omega.col(j), y_j, prior_Pi_AR1[j], j));
    }


    if (verbose) {
      p.increment();
    }
  }

}

