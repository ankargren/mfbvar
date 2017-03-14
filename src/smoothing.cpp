// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#define _USE_MATH_DEFINES // for C++
#include <cmath>

//' @title Smooth and sample from the smoothed distribution
//'
//' @description Functions for smoothing and sampling from the (smoothed) distribution \eqn{p(Z_{1:T}|Y_{1:T}, \Theta)}.
//' @details Implemented in C++.
//' @aliases smoother simulation_smoother generate_mhh loglike
//' @describeIn smoother Compute smoothed states
//' @templateVar Y TRUE
//' @templateVar Lambda TRUE
//' @templateVar Pi_comp TRUE
//' @templateVar Q_comp TRUE
//' @templateVar n_T TRUE
//' @templateVar n_vars TRUE
//' @templateVar n_comp TRUE
//' @templateVar z0 TRUE
//' @templateVar P0 TRUE
//' @template man_template
//' @return For \code{smoother}:
//' \item{}{The smoothed states.}
// [[Rcpp::export]]
arma::mat smoother(           arma::mat Y, arma::mat Lambda, arma::mat Pi_comp, arma::mat Q_comp, int n_T, int n_vars, int n_comp, arma::mat z0, arma::mat P0) {
  /* This function computes the smoothed state vector */
  /****************************************************/
  /* Initialize matrices and cubes */
  arma::mat QQ = Q_comp * Q_comp.t();
  arma::mat mv(n_T, n_vars);
  mv.fill(NA_REAL);
  arma::mat me(n_T, n_vars);
  me.fill(NA_REAL);
  arma::mat mr(n_T, n_comp);
  mr.fill(0);
  arma::mat mu(n_T, n_comp);
  mu.fill(0);
  cube IS(n_vars, n_vars, n_T);
  IS.fill(NA_REAL);
  cube aK(n_comp, n_vars, n_T);
  aK.fill(NA_REAL);
  arma::mat identity_mat(n_comp, n_comp, fill::eye);
  arma::mat YY(n_T, n_vars);
  YY.fill(NA_REAL);
  arma::mat mhh(n_T, n_comp);
  mhh.fill(NA_REAL);
  arma::mat mz = Y.row(0);
  arma::uvec obs_vars = find_finite(mz);

  /* Fill some temporary variables */
  arma::mat h1 = Pi_comp * z0;
  arma::mat P1 = Pi_comp * P0 * Pi_comp.t() + QQ;
  arma::mat mH = Lambda.rows(obs_vars);
  arma::mat vz = mz.cols(obs_vars);

  arma::mat vv = mv.row(0);
  vv.cols(obs_vars) = vz - trans(mH * h1);
  mv.row(0) = vv;

  arma::mat aS = mH * P1 * mH.t();
  arma::mat mIS = IS.slice(0);
  mIS(obs_vars, obs_vars) = inv_sympd(aS);
  IS.slice(0) = mIS;

  arma::mat mK = aK.slice(0);
  mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
  aK.slice(0) = mK;

  arma::mat h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
  arma::mat P2 = (identity_mat - mK.cols(obs_vars) * mH) * P1;

  /* Filtering */
  for (int i = 1; i < n_T; i++) {
    mz = Y.row(i);
    obs_vars = find_finite(mz);

    h1 = Pi_comp * h2;
    P1 = Pi_comp * P2 * Pi_comp.t() + QQ;

    mH = Lambda.rows(obs_vars);
    vz = mz.cols(obs_vars);

    vv = mv.row(i);
    vv.cols(obs_vars) = vz - trans(mH * h1);
    mv.row(i) = vv;

    aS = mH * P1 * mH.t();
    mIS = IS.slice(i);
    mIS(obs_vars, obs_vars) = inv_sympd(aS);
    IS.slice(i) = mIS;

    mK = aK.slice(i);
    mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
    aK.slice(i) = mK;

    h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
    P2 = (identity_mat - mK.cols(obs_vars) * mH) * P1;

  }

  arma::mat ve;
  arma::mat vr;
  arma::mat fk;

  /* Smoothing */
  for (int i = n_T - 1; i >= 1; i--) {
    mz = Y.row(i);
    obs_vars = find_finite(mz);
    mH = Lambda.rows(obs_vars);

    mIS = IS.slice(i);
    mK = aK.slice(i);
    vv = mv.row(i);
    fk = Pi_comp * mK.cols(obs_vars);

    ve = me.row(i);
    ve.cols(obs_vars) = trans(mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)) - fk.t() * trans(mr.row(i)));
    me.row(i) = ve;

    mr.row(i-1) = trans(mH.t() * trans(ve.cols(obs_vars)) + Pi_comp.t() * trans(mr.row(i)));

    mu.row(i) = trans(Q_comp.t() * trans(mr.row(i)));

  }

  mz = Y.row(0);
  obs_vars = find_finite(mz);
  mH = Lambda.rows(obs_vars);

  mIS = IS.slice(0);
  mK = aK.slice(0);
  vv = mv.row(0);
  fk = Pi_comp * mK.cols(obs_vars);

  ve = me.row(0);
  ve.cols(obs_vars) = trans(mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)) - fk.t() * trans(mr.row(0)));
  me.row(0) = ve;

  arma::mat r0 = trans(mH.t() * trans(ve.cols(obs_vars)) + Pi_comp.t() * trans(mr.row(0)));
  mu.row(0) = trans(Q_comp.t() * trans(mr.row(0)));
  arma::mat mu0 = Q_comp.t() * trans(r0);
  mu.insert_rows(0, trans(mu0));

  /* Now get the smoothed states */
  mz = Y.row(0);
  obs_vars = find_finite(mz);
  mH = Lambda.rows(obs_vars);

  arma::mat hz0;
  if (det(P0) == 0) {
    hz0 = z0;
  } else {
    hz0 = z0 + trans(chol(P0)) * as<arma::vec>(rnorm(n_comp));
  }

  mhh.row(0) = trans(Pi_comp * hz0 + Q_comp * trans(mu.row(0)));
  for (int i = 1; i < n_T; i++) {
    mhh.row(i) = trans(Pi_comp * trans(mhh.row(i-1)) + Q_comp * trans(mu.row(i)));
  }

  /* The return is the smoothed state vector */
  return(mhh);
}


//' @describeIn smoother Generate pseudo-state vector
//' @return For \code{generate_mhh}:
//' \item{}{Generated (pseudo-)state vector.}
// [[Rcpp::export]]
arma::mat generate_mhh(       arma::mat Y, arma::mat Lambda, arma::mat Pi_comp, arma::mat Q_comp, int n_T, int n_vars, int n_comp, arma::mat z0, arma::mat P0) {
  /* This function generates the pseudo-state */
  /****************************************************/
  arma::mat mE(n_T, n_comp);
  for (int i = 0; i < n_comp; i++) {
    mE.col(i) = as<arma::vec>(rnorm(n_T));
  }

  arma::mat YY(n_T, n_vars);
  YY.fill(NA_REAL);
  arma::mat mhh(n_T, n_comp);
  mhh.fill(NA_REAL);

  arma::mat hz0;
  if (det(P0) == 0) {
    hz0 = z0;
  } else {
    hz0 = z0 + trans(chol(P0)) * as<arma::vec>(rnorm(n_comp));
  }
  mhh.row(0) = trans(Pi_comp * hz0 + Q_comp * trans(mE.row(0)));
  for (int i = 1; i < n_T; i++) {
    mhh.row(i) = trans(Pi_comp * trans(mhh.row(i-1)) + Q_comp * trans(mE.row(i)));
  }

  return(mhh);
}

//' @describeIn smoother Simulation smoother
//' @return For \code{simulation_smoother}:
//' \item{}{The draw from the posterior distribution.}
// [[Rcpp::export]]
arma::mat simulation_smoother(arma::mat Y, arma::mat Lambda, arma::mat Pi_comp, arma::mat Q_comp, int n_T, int n_vars, int n_comp, arma::mat z0, arma::mat P0) {
  /* This function produces a draw from the posterior distribution */
  /****************************************************/
  arma::mat mhh = generate_mhh(Y, Lambda, Pi_comp, Q_comp, n_T, n_vars, n_comp, z0, P0);
  arma::mat YY(n_T, n_vars);
  YY.fill(NA_REAL);
  arma::mat mz;
  arma::uvec obs_vars;
  arma::mat mH;
  arma::mat mzz;
  mzz.fill(NA_REAL);
  for (int i = 0; i < n_T; i++) {
    mz = Y.row(i);
    obs_vars = find_finite(mz);
    mH = Lambda.rows(obs_vars);
    mzz = YY.row(i);
    mzz.cols(obs_vars) = trans(mH * trans(mhh.row(i)));
    YY.row(i) = mzz;
  }

  arma::mat Z = smoother(Y-YY,  Lambda, Pi_comp, Q_comp, n_T, n_vars, n_comp, z0-z0, P0);
  return(Z+mhh);
  /*return Rcpp::List::create(Rcpp::Named("Z1")  = Z1,
                            Rcpp::Named("Z2")  = Z2,
                            Rcpp::Named("mhh") = mhh,
                            Rcpp::Named("YY") = YY,
                            Rcpp::Named("res") = Z1-Z2+mhh);*/

}


//' @describeIn smoother Compute log-likelihood
//' @return For \code{loglike}:
//' \item{}{An \code{n_T}-long vector of the log-likelihoods. \code{exp(sum(loglike(...)))} is the likelihood.}
// [[Rcpp::export]]
arma::mat loglike(            arma::mat Y, arma::mat Lambda, arma::mat Pi_comp, arma::mat Q_comp, int n_T, int n_vars, int n_comp, arma::mat z0, arma::mat P0) {
  /* This function computes the smoothed state vector */
  /****************************************************/
  /* Initialize matrices and cubes */
  arma::mat QQ = Q_comp * Q_comp.t();
  arma::mat mv(n_T, n_vars);
  mv.fill(NA_REAL);
  arma::mat me(n_T, n_vars);
  me.fill(NA_REAL);
  arma::mat mr(n_T, n_comp);
  mr.fill(0);
  arma::mat mu(n_T, n_comp);
  mu.fill(0);
  cube IS(n_vars, n_vars, n_T);
  IS.fill(NA_REAL);
  cube aK(n_comp, n_vars, n_T);
  aK.fill(NA_REAL);
  arma::mat identity_mat(n_comp, n_comp, fill::eye);
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
  arma::mat P1 = Pi_comp * P0 * Pi_comp.t() + QQ;
  arma::mat mH = Lambda.rows(obs_vars);
  arma::mat vz = mz.cols(obs_vars);

  arma::mat vv = mv.row(0);
  vv.cols(obs_vars) = vz - trans(mH * h1);
  mv.row(0) = vv;

  arma::mat aS = mH * P1 * mH.t();
  arma::mat mIS = IS.slice(0);
  mIS(obs_vars, obs_vars) = inv_sympd(aS);
  IS.slice(0) = mIS;

  arma::mat mK = aK.slice(0);
  mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
  aK.slice(0) = mK;

  arma::mat h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
  arma::mat P2 = (identity_mat - mK.cols(obs_vars) * mH) * P1;

  double log_det_val;
  double log_det_sign;
  /* Filtering */
  for (int i = 1; i < n_T; i++) {
    mz = Y.row(i);
    obs_vars = find_finite(mz);

    h1 = Pi_comp * h2;
    P1 = Pi_comp * P2 * Pi_comp.t() + QQ;

    mH = Lambda.rows(obs_vars);
    vz = mz.cols(obs_vars);

    vv = mv.row(i);
    vv.cols(obs_vars) = vz - trans(mH * h1);
    mv.row(i) = vv;

    aS = mH * P1 * mH.t();
    mIS = IS.slice(i);
    mIS(obs_vars, obs_vars) = inv_sympd(aS);
    IS.slice(i) = mIS;

    mK = aK.slice(i);
    mK.cols(obs_vars) = P1 * mH.t() * mIS(obs_vars, obs_vars);
    aK.slice(i) = mK;

    h2 = h1 + mK.cols(obs_vars) * trans(vv.cols(obs_vars));
    P2 = (identity_mat - mK.cols(obs_vars) * mH) * P1;
    log_det(log_det_val, log_det_sign, aS);
    logl.row(i) = -0.5* obs_vars.n_elem * log(2*M_PI) - (log_det_val + vv.cols(obs_vars) * mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)))*0.5;
  }

  /* The return is the smoothed state vector */
  return(logl);
}
