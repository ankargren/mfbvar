#include "mfbvar.h"
#include "ss_utils.h"
// [[Rcpp::export]]
void mcmc_ss_diffuse(const arma::mat & y_in_p,
                arma::cube& Pi, arma::cube& Sigma, arma::mat& psi, arma::cube& Z,
                arma::cube& Z_fcst, const arma::mat& Lambda_comp, const arma::mat& prior_Pi_Omega,
                const arma::mat& Omega_Pi,
                const arma::mat & D_mat, const arma::mat & dt, const arma::mat & d1,
                const arma::mat & d_fcst_lags, const arma::mat& inv_prior_psi_Omega, const arma::mat& inv_prior_psi_Omega_mean,
                bool check_roots, const arma::mat& Z_1, arma::uword n_reps,
                arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                arma::uword n_T, arma::uword n_fcst, arma::uword n_determ,
                arma::uword n_thin, bool verbose) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps, verbose);

  arma::vec Pi_vec = arma::vec(Pi.begin(), n_vars*(n_vars*n_lags));
  arma::mat Pi_i = Pi.slice(0);
  arma::mat Sigma_i = Sigma.slice(0);
  arma::vec psi_i  = psi.row(0).t();
  arma::mat y_i = y_in_p;
  arma::mat X, post_Pi_Omega_inv, L, b, u1, u2, u4, resid, x;
  arma::mat post_S, mu_mat, mZ, mZ1, mX, Sigma_chol, Sigma_inv;
  arma::mat u3 = arma::vec(n_vars*(n_vars*n_lags));
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

  Sigma_chol = arma::chol(Sigma_i, "lower");
  Sigma_inv = arma::inv_sympd(Sigma_i);
  arma::vec prior_Pi_Omega_vec_inv = 1.0 / arma::vectorise(prior_Pi_Omega);

  // if single freq, we don't need to update
  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_in_p;
  }

  for (arma::uword i = 0; i < n_reps; ++i) {

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
    // Pi
    post_Pi_Omega_inv = arma::kron(Sigma_inv, mX.t() * mX);
    post_Pi_Omega_inv.diag() += prior_Pi_Omega_vec_inv;
    L = arma::chol(post_Pi_Omega_inv, "lower");
    b = arma::vectorise(mX.t() * mZ * Sigma_inv + Omega_Pi);
    u1 = arma::solve(arma::trimatl(L), b);
    u2 = arma::solve(arma::trimatu(L.t()), u1);

    bool stationarity_check = false;
    int num_try = 0, iter = 0;
    double root = 1000;
    while (stationarity_check == false) {
      iter += 1;
      u3.imbue(norm_rand);
      u4 = arma::solve(arma::trimatu(L.t()), u3);
      Pi_vec = u2 + u4;
      Pi_i = arma::trans(arma::reshape(Pi_vec, n_vars*n_lags, n_vars));
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

    resid = mZ - mX * Pi_i.t();
    // Sigma
    post_S = resid.t() * resid;
    Sigma_i = rinvwish(n_T, post_S);
    Sigma_chol = arma::chol(Sigma_i, "lower");
    Sigma_inv = arma::inv_sympd(Sigma_i);

    X = create_X_noint(Z_i, n_lags);
    posterior_psi_iw(psi_i, mu_mat, Pi_i, D_mat, Sigma_i, inv_prior_psi_Omega, mZ + mu_mat, X, inv_prior_psi_Omega_mean, dt, n_determ, n_vars, n_lags);

    arma::vec errors = arma::vec(n_vars);
    if (i % n_thin == 0) {
      if (n_fcst > 0) {
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t() - mu_mat.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors.imbue(norm_rand);
          x = create_X_t_noint(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + Sigma_chol * errors;
        }
        Rcpp::Rcout << "Z_fcst_i: " << arma::size(Z_fcst_i) << std::endl;
        Rcpp::Rcout << "d_fcst_lags: " << arma::size(d_fcst_lags) << std::endl;
        Z_fcst.slice(i/n_thin) = Z_fcst_i.t() + d_fcst_lags * Psi_i.t();
      }
      Z.slice(i/n_thin) = Z_i;
      Sigma.slice(i/n_thin) = Sigma_i;
      Pi.slice(i/n_thin) = Pi_i;
      psi.row(i/n_thin) = psi_i.t();
    }
    if (verbose) {
      p.increment();
    }
  }

}
