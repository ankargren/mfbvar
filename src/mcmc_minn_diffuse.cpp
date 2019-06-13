#include "mfbvar.h"
#include "minn_utils.h"
// [[Rcpp::export]]
void mcmc_minn_diffuse(const arma::mat & y_in_p,
                  arma::cube& Pi, arma::cube& Sigma, arma::cube& Z, arma::cube& Z_fcst,
                  const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                  const arma::mat& Omega_Pi,
                  const arma::mat& Z_1,
                  arma::uword n_reps,
                  arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                  arma::uword n_T, arma::uword n_fcst,
                  arma::uword n_thin, bool verbose) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }


  Progress p(n_reps, verbose);
  arma::vec Pi_vec = arma::vec(Pi.begin(), n_vars*(n_vars*n_lags+1));
  arma::mat Pi_i = Pi.slice(0); //arma::mat(Pi_vec.begin(), n_vars, n_vars*n_lags + 1, false, true);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::mat y_i = y_in_p;
  arma::vec errors = arma::vec(n_vars);
  arma::mat X, post_Pi_Omega_inv, L, b, u1, u2, u4, resid, x;
  arma::mat u3 = arma::vec(n_vars*(n_vars*n_lags + 1));
  arma::mat post_S, Sigma_chol, Sigma_inv;
  arma::mat Z_i = arma::mat(n_lags + y_in_p.n_rows, n_vars, arma::fill::zeros);
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  Z_i.rows(0, n_lags - 1) = Z_1;

  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags);
  }

  Sigma_chol = arma::chol(Sigma_i, "lower");
  Sigma_inv = arma::inv_sympd(Sigma_i);
  arma::vec prior_Pi_Omega_vec_inv = 1.0 / arma::vectorise(prior_Pi_Omega);

  for (arma::uword i = 0; i < n_reps; ++i) {
    if (!single_freq) {
      y_i = simsm_adaptive_cv(y_in_p, Pi_i, Sigma_chol, Lambda_comp, Z_1, n_q, T_b);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags);
    }

    // Pi
    post_Pi_Omega_inv = arma::kron(Sigma_inv, X.t() * X);
    post_Pi_Omega_inv.diag() += prior_Pi_Omega_vec_inv;
    L = arma::chol(post_Pi_Omega_inv, "lower");
    b = arma::vectorise(X.t() * y_i * Sigma_inv + Omega_Pi);
    u1 = arma::solve(arma::trimatl(L), b);
    u2 = arma::solve(arma::trimatu(L.t()), u1);
    u3.imbue(norm_rand);
    u4 = arma::solve(arma::trimatu(L.t()), u3);
    Pi_vec = u2 + u4;
    Pi_i = arma::trans(arma::reshape(Pi_vec, n_vars*n_lags+1, n_vars));
    resid = y_i - X * Pi_i.t(); // Pi_vec and Pi_i use the same memory
    // Sigma
    post_S = resid.t() * resid;
    Sigma_i = rinvwish(n_T, post_S);
    Sigma_chol = arma::chol(Sigma_i, "lower");
    Sigma_inv = arma::inv_sympd(Sigma_i);

    if ((i+1) % n_thin == 0) {
      if (n_fcst > 0) {

        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + Sigma_chol * errors;
        }
        Z_fcst.slice(i/n_thin) = Z_fcst_i.t();
      }

      Z.slice(i/n_thin) = Z_i;
      Sigma.slice(i/n_thin) = Sigma_i;
      Pi.slice(i/n_thin) = Pi_i;
    }
    if (verbose) {
      p.increment();
    }
  }

}


