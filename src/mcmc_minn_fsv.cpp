#include "mfbvar.h"
#include "minn_utils.h"
#include "ss_utils.h"
#include "update_fsv.h"
#include "mvn.h"
#include "mvn_par.h"
#include "update_ng.h"
#include "update_dl.h"
#include <ctime>
// [[Rcpp::export]]
void mcmc_minn_fsv(const arma::mat & y_in_p,
                   arma::cube& Pi, arma::cube& Z, arma::cube& Z_fcst,
                   arma::mat& mu, arma::mat& phi, arma::mat& sigma,
                   arma::cube& f, arma::cube& facload, arma::cube& h,
                   arma::mat & aux, arma::vec & global, arma::mat & local,
                   const arma::mat& Lambda_comp, arma::mat prior_Pi_Omega,
                   const arma::vec& prior_Pi_AR1, const arma::mat& Z_1,
                   double bmu, double Bmu, double a0idi, double b0idi, double a0fac, double b0fac,
                   const Rcpp::NumericVector & Bsigma, double B011inv, double B022inv,
                   const Rcpp::NumericVector & priorh0, const arma::imat & armarestr,
                   const arma::mat & armatau2, // armatau2 is the matrix with prior variance of factor loadings
                   arma::uword n_fac, arma::uword n_reps,
                   arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                   arma::uword n_T, arma::uword n_fcst, arma::uword n_thin, bool verbose,
                   const double a) {
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

  arma::mat eps;
  bool rue = true;
  if ((n_vars*n_lags) > 1.1 * n_T & arma::range(prior_Pi_AR1) < 1e-12) {
    rue = false;
    eps = arma::mat(n_T+n_vars*n_lags+1, n_vars);
  } else {
    eps = arma::mat(n_vars*n_lags+1, n_vars);
  }

  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
    X = create_X(Z_i, n_lags);
  }

  // DL
  bool dl = false;
  double global_i;
  if (a > 0) {
    dl = true;
    global_i = global(0);
  }
  arma::vec aux_i = aux.row(0).t();
  arma::vec local_i = local.row(0).t();

  for (arma::uword i = 0; i < n_reps; ++i) {
    if (!single_freq) {
      Sig_i = arma::exp(0.5 * armah.head_cols(n_vars));
      y_i = simsm_adaptive_univariate(y_in_p, Pi_i, Sig_i, Lambda_comp, Z_1, n_q, T_b, cc_i);
      Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
      X = create_X(Z_i, n_lags);
    }

    y_hat = y_i - X * Pi_i.t();

    if ((i+1) % n_thin == 0) {
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

      global(i/n_thin) = global_i;
      aux.row(i/n_thin) = aux_i.t();
      local.row(i/n_thin) = local_i.t();
    }

    update_fsv(armafacload, armaf, armah, armah0, curpara, armatau2, y_hat.t(),
               bmu, Bmu, a0idi, b0idi, a0fac, b0fac, Bsigma, B011inv, B022inv,
               priorh0, armarestr);

    cc_i = armaf.t() * armafacload.t(); // Common component
    latent_nofac = y_i - cc_i;

    eps.imbue(norm_rand);
    arma::mat output(n_vars*n_lags+1, n_vars);
    if (rue) {
      Pi_parallel_rue Pi_parallel_i(output, latent_nofac, X, prior_Pi_Omega, eps,
                                armah, prior_Pi_AR1, n_T, n_vars, n_lags);
      RcppParallel::parallelFor(0, n_vars, Pi_parallel_i);
    } else {
      Pi_parallel_bcm Pi_parallel_i(output, latent_nofac, X, prior_Pi_Omega, eps,
                                    armah, n_T, n_vars, n_lags);
      RcppParallel::parallelFor(0, n_vars, Pi_parallel_i);
    }

    Pi_i = output.t();

    if (dl) {
      update_dl(prior_Pi_Omega, aux_i, local_i, global_i, Pi_i.t(), n_vars, n_lags, a);
    }


    if (verbose) {
      p.increment();
    }

    if (i % 100 == 0) {
      std::time_t t = std::time(0);   // get time now
      std::tm* now = std::localtime(&t);
      Rcpp::Rcout <<  "Iteration " << i << " at " << now->tm_hour << ':' << now->tm_min <<':'<< now->tm_sec << std::endl;
    }
  }

}

// [[Rcpp::export]]
void mcmc_ss_fsv(const arma::mat & y_in_p,
                 arma::cube& Pi, arma::mat& psi, arma::vec& phi_mu,
                 arma::vec& lambda_mu, arma::mat& omega, arma::cube& Z, arma::cube& Z_fcst,
                 arma::mat& mu, arma::mat& phi, arma::mat& sigma,
                 arma::cube& f, arma::cube& facload, arma::cube& h,
                 const arma::mat& Lambda_comp, arma::mat prior_Pi_Omega,
                 const arma::vec& prior_Pi_AR1,
                 const arma::mat & D_mat, const arma::mat & dt, const arma::mat & d1,
                 const arma::mat & d_fcst_lags, const arma::vec& prior_psi_mean,
                 double c0, double c1, double s, bool check_roots,
                 const arma::mat& Z_1,
                 double bmu, double Bmu, double a0idi, double b0idi, double a0fac, double b0fac,
                 const Rcpp::NumericVector & Bsigma, double B011inv, double B022inv,
                 const Rcpp::NumericVector & priorh0, const arma::imat & armarestr,
                 const arma::mat & armatau2, // armatau2 is the matrix with prior variance of factor loadings
                 arma::uword n_fac, arma::uword n_reps,
                 arma::uword n_q, arma::uword T_b, arma::uword n_lags, arma::uword n_vars,
                 arma::uword n_T, arma::uword n_fcst, arma::uword n_determ, arma::uword n_thin,
                 bool verbose, bool ssng) {
  bool single_freq;
  if (n_q == 0 || n_q == n_vars) {
    single_freq = true;
  } else {
    single_freq = false;
  }

  Progress p(n_reps, verbose);

  arma::mat Pi_i = Pi.slice(0);
  arma::vec psi_i = psi.row(0).t();
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

  arma::mat eps;
  bool rue = true;
  if ((n_vars*n_lags) > 1.1 * n_T & arma::range(prior_Pi_AR1) < 1e-12) {
    rue = false;
    eps = arma::mat(n_T+n_vars*n_lags, n_vars);
  } else {
    eps = arma::mat(n_vars*n_lags, n_vars);
  }

  arma::mat output(n_vars*n_lags, n_vars);


  // ss

  arma::mat mu_mat, mZ, mZ1, mX;
  arma::mat my = arma::mat(arma::size(y_in_p), arma::fill::zeros);
  arma::mat Z_i_demean = arma::mat(n_lags + y_in_p.n_rows, n_vars, arma::fill::zeros);
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
  arma::mat idivar = arma::mat(armah.begin_col(0), armah.n_rows, n_vars, false);

  // ssng
  arma::uword nm = n_vars*n_determ;
  double lambda_mu_i, phi_mu_i, accept, s_prop, M, batch;
  bool adaptive_mh;
  if (ssng) {
    lambda_mu_i = lambda_mu(0);
    phi_mu_i = phi_mu(0);
    accept = 0.0;
    batch = 1.0;

    adaptive_mh = false;

    if (s < 0) {
      M = std::abs(s);
      s = 1.0;
      adaptive_mh = true;
    }
  }

  arma::vec omega_i = omega.row(0).t();
  arma::mat inv_prior_psi_Omega = arma::diagmat(omega_i);
  arma::vec inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
  arma::running_stat<double> stats;

  arma::vec min_vec(2);
  min_vec(0) = 0.01;

  if (single_freq) {
    Z_i.rows(n_lags, n_T + n_lags - 1) = y_i;
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
      Sig_i = arma::exp(0.5 * armah.head_cols(n_vars));
      mZ = simsm_adaptive_univariate(my, Pi_i0, Sig_i, Lambda_comp, mZ1, n_q, T_b, cc_i);
      Z_i.rows(n_lags, n_T + n_lags - 1) = mZ + mu_mat;
    }
    Z_i_demean.rows(0, n_lags - 1) = mZ1;
    Z_i_demean.rows(n_lags, n_T + n_lags - 1) = mZ;

    mX = create_X_noint(Z_i_demean, n_lags);

    y_hat = mZ - mX * Pi_i.t();

    if (i % n_thin == 0) {

      mu_i = curpara_arma.row(0).t();
      phi_i = curpara_arma.row(1).t();
      sigma_i = curpara_arma.row(2).t(); // sigma, not sigma2
      if (n_fcst > 0) {
        vol_pred = armah.tail_rows(1).t();
        Z_fcst_i.head_cols(n_lags) = Z_i.tail_rows(n_lags).t() - mu_mat.tail_rows(n_lags).t();
        for (arma::uword h = 0; h < n_fcst; ++h) {
          errors_sv.imbue(norm_rand);
          errors_var.imbue(norm_rand);
          vol_pred = mu_i + phi_i % (vol_pred - mu_i) + sigma_i % errors_sv; // Twice because we first need it for the volatility, then for the VAR
          error_pred = arma::exp(0.5 * vol_pred) % errors_var;
          x = create_X_t_noint(Z_fcst_i.cols(0+h, n_lags-1+h).t());
          Z_fcst_i.col(n_lags + h) = Pi_i * x + armafacload * error_pred.tail_rows(n_fac) + error_pred.head_rows(n_vars);
        }
        Z_fcst.slice(i/n_thin) = Z_fcst_i.t() + d_fcst_lags * Psi_i.t();
      }


      Z.slice(i/n_thin) = Z_i;
      Pi.slice(i/n_thin) = Pi_i;
      psi.row(i/n_thin) = psi_i.t();

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
    latent_nofac = mZ - cc_i;

    bool stationarity_check = false;
    int num_try = 0, iter = 0;
    double root = 1000;
    while (stationarity_check == false) {
      iter += 1;
      eps.imbue(norm_rand);
      if (rue) {
        Pi_parallel_rue Pi_parallel_i(output, latent_nofac, mX, prior_Pi_Omega, eps,
                                      armah, prior_Pi_AR1, n_T, n_vars, n_lags);
        RcppParallel::parallelFor(0, n_vars, Pi_parallel_i);
      } else {
        Pi_parallel_bcm Pi_parallel_i(output, latent_nofac, mX, prior_Pi_Omega, eps,
                                      armah, n_T, n_vars, n_lags);
        RcppParallel::parallelFor(0, n_vars, Pi_parallel_i);
      }
      Pi_i = output.t();
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

      inv_prior_psi_Omega = arma::diagmat(1/omega_i);
      inv_prior_psi_Omega_mean = prior_psi_mean / omega_i;
    }
    X = create_X_noint(Z_i, n_lags);
    posterior_psi_fsv(psi_i, mu_mat, Pi_i, D_mat, idivar,
                      inv_prior_psi_Omega, Z_i.rows(n_lags, n_T + n_lags - 1), X,
                      armafacload, armaf, inv_prior_psi_Omega_mean,
                      dt, n_determ, n_vars, n_lags);


    if (verbose) {
      p.increment();
    }
  }

}


