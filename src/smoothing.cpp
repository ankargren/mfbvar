#include "mfbvar.h"

#define _USE_MATH_DEFINES // for C++
#include <cmath>

//' @title Smooth and sample from the smoothed distribution
//'
//' @description Functions for smoothing and sampling from the (smoothed) distribution \eqn{p(Z_{1:T}|Y_{1:T}, \Theta)}.
//' @details Implemented in C++.
//' @aliases smoother simulation_smoother generate_mhh loglike
//' @describeIn smoother Compute smoothed states
//' @param Y The data matrix of size \code{n_T x n_vars} with
//' \code{NA} representing missingness. All high-frequency variables must be
//' placed before low-frequency variables.
//' @param Lambda The aggregation matrix (\code{n_vars x (n_vars * n_lags)})
//' @param Pi_comp Matrix with the dynamic coefficients in companion form
//' @param Q_comp The lower-triangular Cholesky decomposition of the covariance matrix (in companion form)
//' @param n_T Number of time periods
//' @param n_vars Number of variables
//' @param n_comp The length of the companion form vector of data (\code{n_vars * n_lags})
//' @param z0 A matrix of size \code{(n_lags * n_vars) x n_vars} of initial values of the latent variable
//' @param P0 The covariance matrix of the initial state (\code{(n_vars * n_lags) x (n_vars * n_lags)})
//' @keywords internal
//' @noRd
//' @return For \code{loglike}:
//' \item{}{An \code{n_T}-long vector of the log-likelihoods. \code{exp(sum(loglike(...)))} is the likelihood.}
// [[Rcpp::export]]
arma::mat loglike(            arma::mat Y, arma::mat Lambda, arma::mat Pi_comp, arma::mat Q_comp, int n_T, int n_vars, int n_comp, arma::mat z0, arma::mat P0) {
  /* This function computes the smoothed state vector */
  /****************************************************/
  /* Initialize matrices and cubes */
  arma::mat QQ = arma::symmatu(Q_comp * Q_comp.t());
  arma::mat mv(n_T, n_vars);
  mv.fill(NA_REAL);
  arma::mat me(n_T, n_vars);
  me.fill(NA_REAL);
  arma::mat mr(n_T, n_comp);
  mr.fill(0);
  arma::mat mu(n_T, n_comp);
  mu.fill(0);
  arma::cube IS(n_vars, n_vars, n_T);
  IS.fill(NA_REAL);
  arma::cube aK(n_comp, n_vars, n_T);
  aK.fill(NA_REAL);
  arma::mat identity_mat(n_comp, n_comp, arma::fill::eye);
  arma::mat YY(n_T, n_vars);
  YY.fill(NA_REAL);
  arma::mat mhh(n_T, n_comp);
  mhh.fill(NA_REAL);
  arma::mat mz = Y.row(0);
  arma::uvec obs_vars = find_finite(mz);
  arma::mat logl(n_T, 1);
  logl.fill(NA_REAL);

  /* Fill some temporary variables */
  arma::mat h1 = Pi_comp * z0;
  arma::mat P1 = arma::symmatu(Pi_comp * P0 * Pi_comp.t() + QQ);
  arma::mat mH = Lambda.rows(obs_vars);
  arma::mat vz = mz.cols(obs_vars);

  arma::mat vv = mv.row(0);
  vv.cols(obs_vars) = vz - trans(mH * h1);
  mv.row(0) = vv;

  arma::mat aS = arma::symmatu(mH * P1 * mH.t());
  arma::mat mIS = IS.slice(0);
  mIS(obs_vars, obs_vars) = inv_sympd(aS);
  IS.slice(0) = mIS;

  arma::mat mK = aK.slice(0);
  mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
  aK.slice(0) = mK;

  arma::mat h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
  arma::mat P2 = arma::symmatu((identity_mat - mK.cols(obs_vars) * mH) * P1);

  double log_det_val;
  double log_det_sign;
  /* Filtering */
  for (int i = 1; i < n_T; i++) {
    mz = Y.row(i);
    obs_vars = find_finite(mz);

    h1 = Pi_comp * h2;
    P1 = arma::symmatu(Pi_comp * P2 * Pi_comp.t() + QQ);

    mH = Lambda.rows(obs_vars);
    vz = mz.cols(obs_vars);

    vv = mv.row(i);
    vv.cols(obs_vars) = vz - trans(mH * h1);
    mv.row(i) = vv;

    aS = arma::symmatu(mH * P1 * mH.t());
    mIS = IS.slice(i);
    mIS(obs_vars, obs_vars) = inv_sympd(aS);
    IS.slice(i) = mIS;

    mK = aK.slice(i);
    mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
    aK.slice(i) = mK;

    h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
    P2 = arma::symmatu((identity_mat - mK.cols(obs_vars) * mH) * P1);
    log_det(log_det_val, log_det_sign, aS);
    logl.row(i) = -0.5* obs_vars.n_elem * std::log(2*M_PI) - (log_det_val + vv.cols(obs_vars) * mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)))*0.5;
  }

  /* The return is the smoothed state vector */
  return(logl);
}
