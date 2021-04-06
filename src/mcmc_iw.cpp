#include "mfbvar.h"
#include "minn_utils.h"
#include "ss_utils.h"
#include "update_ng.h"
#include "update_iw.h"
// [[Rcpp::export]]
void mcmc_minn_iw(const arma::mat & y_in_p,
                  arma::cube& Pi, arma::cube& Sigma, arma::cube& Z, arma::cube& Z_fcst,
                  const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                  const arma::mat& inv_prior_Pi_Omega,
                  const arma::mat& Omega_Pi, const arma::mat& prior_Pi_mean,
                  const arma::mat & prior_S, bool check_roots, const arma::mat& Z_1,
                  arma::uword n_reps, arma::uword n_burnin,
                  arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                  arma::uword n_T, arma::uword n_fcst,
                  arma::uword n_thin, bool verbose, int prior_nu,
                  bool fixate_Pi, bool fixate_Sigma, bool fixate_Z) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }


  Progress p(n_reps+n_burnin, verbose);
  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::mat y_i = Z.slice(0).rows(n_lags, n_T + n_lags - 1);
  arma::mat Z_i = Z.slice(0);

  arma::vec errors = arma::vec(n_vars);
  arma::mat X, post_Pi_Omega, post_Pi;
  arma::mat post_S, Sigma_chol, x;
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  int post_nu = n_T + prior_nu;

  if (single_freq) {
    y_i = y_in_p;
  }
  if (single_freq || fixate_Z) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags);
    update_iw(post_Pi_Omega, post_Pi, post_S, X, y_i, prior_Pi_mean,
              prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S);
  }

  Sigma_chol = arma::chol(Sigma_i, "lower");

  for (arma::uword i = 0; i < n_reps + n_burnin; ++i) {
    if (!single_freq && !fixate_Z) {
      y_i = simsm_adaptive_cv(y_in_p, Pi_i, Sigma_chol, Lambda_comp, Z_1, n_q, T_b);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags);
      update_iw(post_Pi_Omega, post_Pi, post_S, X, y_i, prior_Pi_mean,
                prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S);
    }
    if (!fixate_Sigma) {
      Sigma_i = rinvwish(post_nu, post_S);
      Sigma_chol = arma::chol(Sigma_i, "lower");
    }
    if (!fixate_Pi) {
      Pi_i = rmatn(post_Pi.t(), post_Pi_Omega, Sigma_i);
    }

    if ((i+1) % n_thin == 0 && i >= n_burnin) {
      if (n_fcst > 0) {
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + Sigma_chol * errors;
        }
        Z_fcst.slice((i-n_burnin)/n_thin) = Z_fcst_i.t();
      }

      Z.slice((i-n_burnin)/n_thin) = Z_i;
      Sigma.slice((i-n_burnin)/n_thin) = Sigma_i;
      Pi.slice((i-n_burnin)/n_thin) = Pi_i;
    }
    if (verbose) {
      p.increment();
    }
  }

}

// [[Rcpp::export]]
void mcmc_ssng_iw(const arma::mat & y_in_p,
                  arma::cube& Pi, arma::cube& Sigma, arma::cube& psi, arma::cube& phi_mu,
                  arma::cube& lambda_mu, arma::cube& omega, arma::cube& Z,
                  arma::cube& Z_fcst, const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                  const arma::mat& inv_prior_Pi_Omega,
                  const arma::mat& Omega_Pi, const arma::mat& prior_Pi_mean,
                  const arma::mat & prior_S,
                  const arma::mat & D_mat, const arma::mat & dt, const arma::mat & d1,
                  const arma::mat & d_fcst_lags, const arma::vec& prior_psi_mean,
                  double c0, double c1, double s,
                  bool check_roots, const arma::mat& Z_1, arma::uword n_reps, arma::uword n_burnin,
                  arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                  arma::uword n_T, arma::uword n_fcst, arma::uword n_determ,
                  arma::uword n_thin, bool verbose, int prior_nu, bool ssng,
                  bool fixate_Pi, bool fixate_Sigma, bool fixate_Z, bool fixate_psi,
                  bool fixate_phi_mu, bool fixate_lambda_mu, bool fixate_omega) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps+n_burnin, verbose);

  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::vec psi_i  = psi.slice(0);
  arma::mat Z_i = Z.slice(0);
  double lambda_mu_i = arma::as_scalar(lambda_mu.slice(0));
  double phi_mu_i = arma::as_scalar(phi_mu.slice(0));
  arma::vec omega_i = omega.slice(0);

  arma::mat X, post_Pi_Omega, post_Pi;
  arma::mat post_S, x, mu_mat, mZ, mZ1, mX;
  arma::mat my = arma::mat(arma::size(y_in_p), arma::fill::zeros);

  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  arma::mat Z_i_demean = Z_i;
  Z_i.rows(0, n_lags - 1) = Z_1;

  arma::mat Pi_i0 = arma::mat(n_vars, n_vars*n_lags+1, arma::fill::zeros);
  int post_nu = n_T + prior_nu;
  arma::mat Sigma_chol = arma::chol(Sigma_i, "lower");

  arma::mat Psi_i = arma::mat(psi_i.begin(), n_vars, n_determ, false, true);
  mu_mat = dt * Psi_i.t();
  arma::uword n_Lambda = Lambda_comp.n_cols/Lambda_comp.n_rows;
  arma::mat mu_long = arma::mat(n_Lambda+n_T, n_vars, arma::fill::zeros);
  arma::rowvec Lambda_single = arma::rowvec(n_Lambda, arma::fill::zeros);
  for (arma::uword i = 0; i < n_Lambda; ++i) {
    Lambda_single(i) = Lambda_comp.at(0, i*n_q);
  }

  arma::uword nm = n_vars*n_determ;
  arma::mat inv_prior_psi_Omega = arma::diagmat(1.0/omega_i);
  arma::vec inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;

  double M;
  arma::running_stat<double> stats;
  double accept = 0.0;
  bool adaptive_mh = false;
  if (s < 0) {
    M = std::abs(s);
    s = 1.0;
    adaptive_mh = true;
  }

  // if single freq, we don't need to update
  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_in_p;
  }
  for (arma::uword i = 0; i < n_reps + n_burnin; ++i) {

    if (!single_freq) {
      update_demean(my, mu_long, y_in_p, mu_mat, d1, Psi_i, Lambda_single, n_vars,
                    n_q, n_Lambda, n_T);
    } else {
      mZ = y_in_p - mu_mat;
    }

    mZ1 = Z_1 - d1 * Psi_i.t();
    Pi_i0.cols(1, n_vars*n_lags) = Pi_i;

    if (!single_freq) {
      if (!fixate_Z) {
        mZ = simsm_adaptive_cv(my, Pi_i0, Sigma_chol, Lambda_comp, mZ1, n_q, T_b);
        Z_i.rows(n_lags, n_T + n_lags - 1) = mZ + mu_mat;
      } else {
        mZ = Z_i.rows(n_lags, n_T + n_lags - 1) - mu_mat;
      }
    }

    Z_i_demean.rows(0, n_lags - 1) = mZ1;
    Z_i_demean.rows(n_lags, n_T + n_lags - 1) = mZ;

    mX = create_X_noint(Z_i_demean, n_lags);
    update_iw(post_Pi_Omega, post_Pi, post_S, mX, mZ, prior_Pi_mean,
              prior_Pi_Omega, inv_prior_Pi_Omega, Omega_Pi, prior_S);
    if (!fixate_Sigma) {
      Sigma_i = rinvwish(post_nu, arma::symmatu(post_S));
      Sigma_chol = arma::chol(Sigma_i, "lower");
    }
    if (!fixate_Pi) {
      Pi_i = sample_Pi_mat(post_Pi, post_Pi_Omega, Sigma_i, check_roots, n_vars);
    }

    if (ssng) {
      update_ng(phi_mu_i, lambda_mu_i, omega_i, nm, c0, c1, s, psi_i,
                prior_psi_mean, accept, fixate_phi_mu, fixate_lambda_mu,
                fixate_omega);
      if (adaptive_mh) {
        update_s(s, stats, accept, i, M);
      }
      inv_prior_psi_Omega = arma::diagmat(1.0/omega_i);
      inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
    }

    X = create_X_noint(Z_i, n_lags);
    if (!fixate_psi) {
      posterior_psi_iw(psi_i, mu_mat, Pi_i, D_mat, Sigma_i, inv_prior_psi_Omega,
                       mZ + mu_mat, X, inv_prior_psi_Omega_mean, dt, n_determ,
                       n_vars, n_lags);
    }

    arma::vec errors = arma::vec(n_vars);
    if ((i+1) % n_thin == 0 && i>=n_burnin) {
      if (n_fcst > 0) {
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t() - mu_mat.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t_noint(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + Sigma_chol * errors;
        }
        Z_fcst.slice((i-n_burnin)/n_thin) = Z_fcst_i.t() + d_fcst_lags * Psi_i.t();
      }
      Z.slice((i-n_burnin)/n_thin) = Z_i;
      Sigma.slice((i-n_burnin)/n_thin) = Sigma_i;
      Pi.slice((i-n_burnin)/n_thin) = Pi_i;
      psi.slice((i-n_burnin)/n_thin) = psi_i;
      if (ssng) {
        phi_mu.slice((i-n_burnin)/n_thin) = phi_mu_i;
        lambda_mu.slice((i-n_burnin)/n_thin) = lambda_mu_i;
        omega.slice((i-n_burnin)/n_thin) = omega_i;
      }

    }
    if (verbose) {
      p.increment();
    }
  }

}

