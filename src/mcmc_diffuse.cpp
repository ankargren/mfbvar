#include "mfbvar.h"
#include "minn_utils.h"
#include "ss_utils.h"
#include "update_ng.h"
#include "update_dl.h"
// [[Rcpp::export]]
void mcmc_minn_diffuse(const arma::mat & y_in_p,
                  arma::cube& Pi, arma::cube& Sigma, arma::cube& Z, arma::cube& Z_fcst,
                  arma::cube& aux, arma::cube& global, arma::cube& local,
                  arma::cube& slice,
                  const arma::mat& Lambda, arma::mat prior_Pi_Omega,
                  arma::vec prior_Pi_mean_vec, bool check_roots,
                  const arma::mat& Z_1,
                  arma::uword n_reps, arma::uword n_burnin,
                  arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                  arma::uword n_T, arma::uword n_fcst,
                  arma::uword n_thin, bool verbose,
                  const double a, bool gig,
                  bool fixate_Pi, bool fixate_Sigma, bool fixate_Z,
                  bool fixate_aux, bool fixate_global, bool fixate_local) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps + n_burnin, verbose);
  arma::vec Pi_vec = arma::vec(Pi.begin(), n_vars*(n_vars*n_lags+1));
  arma::mat Pi_i = Pi.slice(0); //arma::mat(Pi_vec.begin(), n_vars, n_vars*n_lags + 1, false, true);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::mat Z_i = Z.slice(0);
  arma::vec errors = arma::vec(n_vars);
  arma::mat X, post_Pi_Omega_inv, L, b, u1, u2, u4, resid, x, y_i;
  arma::mat u3 = arma::vec(n_vars*(n_vars*n_lags + 1));
  arma::mat post_S, Sigma_chol, Sigma_inv;
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  Z_i.rows(0, n_lags - 1) = Z_1;

  if (single_freq) {
    y_i = y_in_p;
  } else {
    y_i = Z.slice(0).rows(n_lags, n_T + n_lags - 1);
  }

  if (single_freq || fixate_Z) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags, true);
  }

  // DL
  bool dl = false;
  double global_i;
  if (a > 0) {
    dl = true;
    global_i = arma::as_scalar(global.slice(0));
  }
  arma::vec aux_i = aux.slice(0);
  arma::vec local_i = local.slice(0);
  arma::vec slice_i = slice.slice(0);

  if (dl) {
    prior_Pi_Omega.rows(1, n_vars*n_lags) =
      arma::reshape(aux_i % arma::pow(global_i * local_i, 2.0),
                    n_vars*n_lags, n_vars);
  }

  Sigma_chol = arma::chol(Sigma_i, "lower");
  Sigma_inv = arma::inv_sympd(Sigma_i);
  arma::vec prior_Pi_Omega_vec_inv = 1.0 / arma::vectorise(prior_Pi_Omega);
  arma::vec Omega_Pi = prior_Pi_mean_vec % prior_Pi_Omega_vec_inv;

  for (arma::uword i = 0; i < n_reps + n_burnin; ++i) {
    if (!single_freq && !fixate_Z) {
      y_i = simsm_adaptive_cv(y_in_p, Pi_i, Sigma_chol,
                              Lambda, Z_1, n_q, T_b);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags, true);
    }

    // Pi
    if (!fixate_Pi) {
      Pi_i = sample_Pi_vec(Sigma_inv, X, y_i, prior_Pi_Omega_vec_inv, Omega_Pi,
                           check_roots, n_vars, n_lags);
    }
    resid = y_i - X * Pi_i.t();

    // Sigma
    if (!fixate_Sigma) {
      post_S = resid.t() * resid;
      Sigma_i = rinvwish(n_T, post_S);
      Sigma_chol = arma::chol(Sigma_i, "lower");
      Sigma_inv = arma::inv_sympd(Sigma_i);
    }

    if (dl) {
      update_dl(prior_Pi_Omega, aux_i, local_i, global_i, Pi_i.t(), n_vars,
                n_lags, a, slice_i, fixate_aux, fixate_global, fixate_local,
                gig, true);
      prior_Pi_Omega_vec_inv = 1.0 / arma::vectorise(prior_Pi_Omega);
      Omega_Pi = prior_Pi_mean_vec % prior_Pi_Omega_vec_inv;
    }

    if (((i+1) % n_thin == 0) && (i >= n_burnin)) {
      if (n_fcst > 0) {

        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t(Z_fcst_i.cols(0+h, n_lags-1+h).t(), true);
          Z_fcst_i.col(n_lags + h) = Pi_i * x + Sigma_chol * errors;
        }
        Z_fcst.slice((i-n_burnin)/n_thin) = Z_fcst_i.t();
      }

      Z.slice((i-n_burnin)/n_thin) = Z_i;
      Sigma.slice((i-n_burnin)/n_thin) = Sigma_i;
      Pi.slice((i-n_burnin)/n_thin) = Pi_i;

      if (dl) {
        global.slice((i-n_burnin)/n_thin) = global_i;
        aux.slice((i-n_burnin)/n_thin) = aux_i;
        local.slice((i-n_burnin)/n_thin) = local_i;
        if (!gig) {
          slice.slice((i-n_burnin)/n_thin) = slice_i;
        }
      }
    }
    if (verbose) {
      p.increment();
    }
  }

}



// [[Rcpp::export]]
void mcmc_ssng_diffuse(const arma::mat & y_in_p,
                       arma::cube& Pi, arma::cube& Sigma, arma::cube& psi,
                       arma::cube& phi_mu, arma::cube& lambda_mu,
                       arma::cube& omega, arma::cube& Z, arma::cube& Z_fcst,
                       const arma::mat& Lambda,
                       const arma::mat& prior_Pi_Omega,
                       const arma::vec & prior_Pi_mean_vec,
                       const arma::mat & D_mat, const arma::mat & dt,
                       const arma::mat & d1, const arma::mat & d_fcst_lags,
                       const arma::vec& prior_psi_mean, double c0, double c1,
                       double s, bool check_roots, const arma::mat& Z_1,
                       arma::uword n_reps, arma::uword n_burnin,
                       arma::uword n_q, arma::uword T_b, arma::uword n_lags,
                       arma::uword n_vars, arma::uword n_T, arma::uword n_fcst,
                       arma::uword n_determ, arma::uword n_thin, bool verbose,
                       bool ssng, bool fixate_Pi, bool fixate_Sigma,
                       bool fixate_Z, bool fixate_psi, bool fixate_phi_mu,
                       bool fixate_lambda_mu, bool fixate_omega) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps + n_burnin, verbose);
  arma::vec Pi_vec = arma::vec(Pi.begin(), n_vars*(n_vars*n_lags));
  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::vec psi_i  = psi.slice(0);
  arma::mat y_i = Z.slice(0).rows(n_lags, n_T + n_lags - 1);
  arma::mat Z_i = Z.slice(0);
  double lambda_mu_i = arma::as_scalar(lambda_mu.slice(0));
  double phi_mu_i = arma::as_scalar(phi_mu.slice(0));
  arma::vec omega_i = omega.slice(0);

  arma::mat X, resid, x;
  arma::mat post_S, mu_mat, mZ, mZ1, mX, Sigma_chol, Sigma_inv;
  arma::mat my = arma::mat(arma::size(y_in_p), arma::fill::zeros);

  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  arma::mat Z_i_demean = Z_i;
  Z_i.rows(0, n_lags - 1) = Z_1;

  arma::mat Pi_i0 = arma::mat(n_vars, n_vars*n_lags+1, arma::fill::zeros);
  arma::mat Pi_comp = arma::mat(n_vars*n_lags, n_vars*n_lags, arma::fill::zeros);
  Pi_comp.submat(n_vars, 0, n_vars*n_lags - 1, n_vars*(n_lags-1) - 1) = arma::eye(n_vars*(n_lags-1), n_vars*(n_lags-1));

  arma::mat Psi_i = arma::mat(psi_i.begin(), n_vars, n_determ, false, true);
  mu_mat = dt * Psi_i.t();
  arma::uword n_Lambda = Lambda.n_cols/Lambda.n_rows;
  arma::mat mu_long = arma::mat(n_Lambda+n_T, n_vars, arma::fill::zeros);
  arma::rowvec Lambda_single = arma::rowvec(n_Lambda, arma::fill::zeros);
  for (arma::uword i = 0; i < n_Lambda; ++i) {
    Lambda_single(i) = Lambda.at(0, i*n_q);
  }

  Sigma_chol = arma::chol(Sigma_i, "lower");
  Sigma_inv = arma::inv_sympd(Sigma_i);
  arma::vec prior_Pi_Omega_vec_inv = 1.0 / arma::vectorise(prior_Pi_Omega);
  arma::vec Omega_Pi = prior_Pi_mean_vec % prior_Pi_Omega_vec_inv;

  // if single freq, we don't need to update
  if (single_freq) {
    y_i = y_in_p;
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_in_p;
  }

  // NG stuff
  arma::uword nm = n_vars*n_determ;
  arma::mat inv_prior_psi_Omega = arma::diagmat(1.0/omega_i);
  arma::vec inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
  double M;
  arma::running_stat<double> stats;
  arma::uword accept = 0;
  bool adaptive_mh = false;
  if (s < 0) {
    M = std::abs(s);
    s = 1.0;
    adaptive_mh = true;
  }
  for (arma::uword i = 0; i < n_reps + n_burnin; ++i) {
    if (!single_freq) {
      update_demean(my, mu_long, y_in_p, mu_mat, d1, Psi_i, Lambda_single, n_vars,
                    n_q, n_Lambda, n_T);
    } else {
      // Even if single freq, mZ needs to be updated
      mZ = y_in_p - mu_mat;
    }

    mZ1 = Z_1 - d1 * Psi_i.t();
    Pi_i0.cols(1, n_vars*n_lags) = Pi_i;

    if (!single_freq) {
      if (!fixate_Z) {
        mZ = simsm_adaptive_cv(my, Pi_i0, Sigma_chol, Lambda, mZ1, n_q, T_b);
        Z_i.rows(n_lags, n_T + n_lags - 1) = mZ + mu_mat;
      } else {
        mZ = Z_i.rows(n_lags, n_T + n_lags - 1) - mu_mat;
      }
    }

    Z_i_demean.rows(0, n_lags - 1) = mZ1;
    Z_i_demean.rows(n_lags, n_T + n_lags - 1) = mZ;

    mX = create_X(Z_i_demean, n_lags, false);
    if (!fixate_Pi) {
      Pi_i = sample_Pi_vec(Sigma_inv, mX, mZ, prior_Pi_Omega_vec_inv, Omega_Pi,
                           check_roots, n_vars, n_lags);
    }
    resid = mZ - mX * Pi_i.t();

    // Sigma
    if (!fixate_Sigma) {
      post_S = resid.t() * resid;
      Sigma_i = rinvwish(n_T, post_S);
      Sigma_chol = arma::chol(Sigma_i, "lower");
      Sigma_inv = arma::inv_sympd(Sigma_i);
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

    X = create_X(Z_i, n_lags, false);
    if (!fixate_psi) {
      posterior_psi_iw(psi_i, mu_mat, Pi_i, D_mat, Sigma_i, inv_prior_psi_Omega,
                       mZ + mu_mat, X, inv_prior_psi_Omega_mean, dt, n_determ,
                       n_vars, n_lags);
    }

    arma::vec errors = arma::vec(n_vars);
    if (((i+1) % n_thin == 0) && (i >= n_burnin)) {
      if (n_fcst > 0) {
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t() - mu_mat.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t(Z_fcst_i.cols(0+h, n_lags-1+h).t(), false);
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


