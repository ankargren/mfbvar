#include "mfbvar.h"
#include "minn_utils.h"
#include "ss_utils.h"
#include "update_ng.h"
// [[Rcpp::export]]
void mcmc_minn_iw(const arma::mat & y_in_p,
                  arma::cube& Pi, arma::cube& Sigma, arma::cube& Z, arma::cube& Z_fcst,
                  const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                  const arma::mat& inv_prior_Pi_Omega,
                  const arma::mat& Omega_Pi, const arma::mat& prior_Pi_mean,
                  const arma::mat & prior_S, const arma::mat& Z_1,
                  arma::uword n_reps, arma::uword n_burnin,
                  arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                  arma::uword n_T, arma::uword n_fcst,
                  arma::uword n_thin, bool verbose, int prior_nu) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }


  Progress p(n_reps+n_burnin, verbose);
  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::mat y_i = y_in_p;
  arma::vec errors = arma::vec(n_vars);
  arma::mat X, XX, XX_inv, Pi_sample, post_Pi_Omega, post_Pi;
  arma::mat S, Pi_diff, post_S, Sigma_chol, x;
  arma::mat Z_i = arma::mat(n_lags + y_in_p.n_rows, n_vars, arma::fill::zeros);
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  Z_i.rows(0, n_lags - 1) = Z_1;
  int post_nu = n_T + n_vars + prior_nu;

  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags);
    XX = X.t() * X;
    XX_inv = arma::inv_sympd(XX);
    Pi_sample = XX_inv * (X.t() * y_i);
    post_Pi_Omega = arma::inv_sympd(inv_prior_Pi_Omega + XX);
    post_Pi = post_Pi_Omega * (Omega_Pi + X.t() * y_i);
    S = arma::trans((y_i - X * Pi_sample)) * (y_i - X * Pi_sample);
    Pi_diff = prior_Pi_mean - Pi_sample;
    post_S = prior_S + S + Pi_diff.t() * arma::inv_sympd(prior_Pi_Omega + XX_inv) * Pi_diff;
  }

  Sigma_chol = arma::chol(Sigma_i, "lower");

  for (arma::uword i = 0; i < n_reps + n_burnin; ++i) {
    if (!single_freq) {
      y_i = simsm_adaptive_cv(y_in_p, Pi_i, Sigma_chol, Lambda_comp, Z_1, n_q, T_b);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags);
      XX = X.t() * X;
      XX_inv = arma::inv_sympd(XX);
      Pi_sample = XX_inv * (X.t() * y_i);
      post_Pi_Omega = arma::inv_sympd(inv_prior_Pi_Omega + XX);
      post_Pi = post_Pi_Omega * (Omega_Pi + X.t() * y_i);
      S = arma::trans((y_i - X * Pi_sample)) * (y_i - X * Pi_sample);
      Pi_diff = prior_Pi_mean - Pi_sample;
      post_S = prior_S + S + Pi_diff.t() * arma::inv_sympd(prior_Pi_Omega + XX_inv) * Pi_diff;
    }
    Sigma_i = rinvwish(post_nu, post_S);
    Sigma_chol = arma::chol(Sigma_i, "lower");
    Pi_i = rmatn(post_Pi.t(), post_Pi_Omega, Sigma_i);

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
                  arma::cube& Pi, arma::cube& Sigma, arma::mat& psi, arma::vec& phi_mu,
                  arma::vec& lambda_mu, arma::mat& omega, arma::cube& Z,
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
                  arma::uword n_thin, bool verbose, bool ssng) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps+n_burnin, verbose);

  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::vec psi_i  = psi.row(0).t();
  arma::mat y_i, X, XX, XX_inv, Pi_sample, post_Pi_Omega, post_Pi;
  arma::mat S, Pi_diff, post_S, x, mu_mat, mZ, mZ1, mX;
  arma::mat my = arma::mat(arma::size(y_in_p), arma::fill::zeros);

  arma::mat Z_i = arma::mat(n_lags + y_in_p.n_rows, n_vars, arma::fill::zeros);
  arma::mat Z_fcst_i = arma::mat(n_vars, n_lags + n_fcst);
  arma::mat Z_i_demean = Z_i;
  Z_i.rows(0, n_lags - 1) = Z_1;

  arma::mat Pi_i0 = arma::mat(n_vars, n_vars*n_lags+1, arma::fill::zeros);
  arma::mat Pi_comp = arma::mat(n_vars*n_lags, n_vars*n_lags, arma::fill::zeros);
  Pi_comp.submat(n_vars, 0, n_vars*n_lags - 1, n_vars*(n_lags-1) - 1) = arma::eye(n_vars*(n_lags-1), n_vars*(n_lags-1));

  arma::mat Psi_i = arma::mat(psi_i.begin(), n_vars, n_determ, false, true);
  mu_mat = dt * Psi_i.t();
  arma::uword n_Lambda = Lambda_comp.n_cols/Lambda_comp.n_rows;
  arma::mat mu_long = arma::mat(n_Lambda+n_T, n_vars, arma::fill::zeros);
  arma::rowvec Lambda_single = arma::rowvec(n_Lambda, arma::fill::zeros);
  for (arma::uword i = 0; i < n_Lambda; ++i) {
    Lambda_single(i) = Lambda_comp.at(0, i*n_q);
  }

  int post_nu = n_T + n_vars + 2;
  arma::mat Sigma_chol = arma::chol(Sigma_i, "lower");

  arma::uword nm = n_vars*n_determ;
  double lambda_mu_i = lambda_mu(0);
  double phi_mu_i = phi_mu(0);
  arma::vec omega_i = omega.row(0).t();
  arma::mat inv_prior_psi_Omega = arma::diagmat(1.0/omega_i);
  arma::vec inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
  double M, batch = 1.0;
  arma::running_stat<double> stats;
  double accept = 0.0;
  bool adaptive_mh = false;
  double s_prop;
  if (s < 0) {
    M = std::abs(s);
    s = 1.0;
    adaptive_mh = true;
  }
  arma::vec min_vec(2);
  min_vec(0) = 0.01;

  // if single freq, we don't need to update
  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_in_p;
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
      mZ = simsm_adaptive_cv(my, Pi_i0, Sigma_chol, Lambda_comp, mZ1, n_q, T_b);
      Z_i.rows(n_lags, n_T + n_lags - 1) = mZ + mu_mat;
    }

    Z_i_demean.rows(0, n_lags - 1) = mZ1;
    Z_i_demean.rows(n_lags, n_T + n_lags - 1) = mZ;

    mX = create_X_noint(Z_i_demean, n_lags);
    XX = mX.t() * mX;
    XX_inv = arma::inv_sympd(XX);
    Pi_sample = XX_inv * (mX.t() * mZ);

    post_Pi_Omega = arma::inv_sympd(inv_prior_Pi_Omega + XX);
    post_Pi = post_Pi_Omega * (Omega_Pi + mX.t() * mZ);
    S = arma::trans((mZ - mX * Pi_sample)) * (mZ - mX * Pi_sample);
    Pi_diff = prior_Pi_mean - Pi_sample;
    post_S = prior_S + S + Pi_diff.t() * arma::inv_sympd(prior_Pi_Omega + XX_inv) * Pi_diff;
    Sigma_i = rinvwish(post_nu, arma::symmatu(post_S)); //Fixed in 9.400.3
    Sigma_chol = arma::chol(Sigma_i, "lower");
    bool stationarity_check = false;
    int num_try = 0, iter = 0;
    double root = 1000;
    while (stationarity_check == false) {
      iter += 1;
      Pi_i = rmatn(post_Pi.t(), post_Pi_Omega, Sigma_i);
      if (check_roots) {
        Pi_comp.rows(0, n_vars-1) = Pi_i;
        root = max_eig_cpp(Pi_comp);
      } else {
        root = 0.0;
      }
      if (root < 1.0) {
        stationarity_check = true;
        num_try = iter;
      }
      if (iter == 1000) {
        Rcpp::stop("Attemped to draw stationary Pi 1,000 times.");
      }
    }

    if (ssng) {
      update_ng(phi_mu_i, lambda_mu_i, omega_i, nm, c0, c1, s, psi_i, prior_psi_mean, accept);
      if (adaptive_mh) {
        stats(accept);
        if (i % 100 == 0) {
          batch += 1.0;
          min_vec(1) = std::pow(batch, -0.5);
          if (stats.mean() > 0.44) {
            s_prop = log(s) + arma::min(min_vec);
            if (s_prop < M){
              s = std::exp(s_prop);
            }
          } else {
            s_prop = log(s) - arma::min(min_vec);
            if (s_prop > -M){
              s = std::exp(s_prop);
            }
          }
          stats.reset();
        }
      }

      inv_prior_psi_Omega = arma::diagmat(1.0/omega_i);
      inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
    }


    X = create_X_noint(Z_i, n_lags);

    posterior_psi_iw(psi_i, mu_mat, Pi_i, D_mat, Sigma_i, inv_prior_psi_Omega, mZ + mu_mat, X, inv_prior_psi_Omega_mean, dt, n_determ, n_vars, n_lags);
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
      psi.row((i-n_burnin)/n_thin) = psi_i.t();
      if (ssng) {
        phi_mu((i-n_burnin)/n_thin) = phi_mu_i;
        lambda_mu((i-n_burnin)/n_thin) = lambda_mu_i;
        omega.row((i-n_burnin)/n_thin) = omega_i.t();
      }

    }
    if (verbose) {
      p.increment();
    }
  }

}

