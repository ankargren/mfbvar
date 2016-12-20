// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat smoother(           arma::mat mZ, arma::mat Lambda, arma::mat mF, arma::mat mQ, int iT, int ip, int iq, arma::mat h0, arma::mat P0) {
  /* This function computes the smoothed state vector */
  /****************************************************/
  /* Initialize matrices and cubes */
  arma::mat QQ = mQ * mQ.t();
  arma::mat mv(iT, ip);
  mv.fill(NA_REAL);
  arma::mat me(iT, ip);
  me.fill(NA_REAL);
  arma::mat mr(iT, iq);
  mr.fill(0);
  arma::mat mu(iT, iq);
  mu.fill(0);
  cube IS(ip, ip, iT);
  IS.fill(NA_REAL);
  cube aK(iq, ip, iT);
  aK.fill(NA_REAL);
  arma::mat identity_mat(iq, iq, fill::eye);
  arma::mat mZZ(iT, ip);
  mZZ.fill(NA_REAL);
  arma::mat mhh(iT, iq);
  mhh.fill(NA_REAL);
  arma::mat mz = mZ.row(0);
  arma::uvec obs_vars = find_finite(mz);

  /* Fill some temporary variables */
  arma::mat h1 = mF * h0;
  arma::mat P1 = mF * P0 * mF.t() + QQ;
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
  for (int i = 1; i < iT; i++) {
   mz = mZ.row(i);
   obs_vars = find_finite(mz);

   h1 = mF * h2;
   P1 = mF * P2 * mF.t() + QQ;

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
  for (int i = iT - 1; i >= 1; i--) {
   mz = mZ.row(i);
   obs_vars = find_finite(mz);
   mH = Lambda.rows(obs_vars);

   mIS = IS.slice(i);
   mK = aK.slice(i);
   vv = mv.row(i);
   fk = mF * mK.cols(obs_vars);

   ve = me.row(i);
   ve.cols(obs_vars) = trans(mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)) - fk.t() * trans(mr.row(i)));
   me.row(i) = ve;

   mr.row(i-1) = trans(mH.t() * trans(ve.cols(obs_vars)) + mF.t() * trans(mr.row(i)));

   mu.row(i) = trans(mQ.t() * trans(mr.row(i)));

  }

  mz = mZ.row(0);
  obs_vars = find_finite(mz);
  mH = Lambda.rows(obs_vars);

  mIS = IS.slice(0);
  mK = aK.slice(0);
  vv = mv.row(0);
  fk = mF * mK.cols(obs_vars);

  ve = me.row(0);
  ve.cols(obs_vars) = trans(mIS(obs_vars, obs_vars) * trans(vv.cols(obs_vars)) - fk.t() * trans(mr.row(0)));
  me.row(0) = ve;

  arma::mat r0 = trans(mH.t() * trans(ve.cols(obs_vars)) + mF.t() * trans(mr.row(0)));
  mu.row(0) = trans(mQ.t() * trans(mr.row(0)));
  arma::mat mu0 = mQ.t() * trans(r0);
  mu.insert_rows(0, trans(mu0));

  /* Now get the smoothed states */
  mz = mZ.row(0);
  obs_vars = find_finite(mz);
  mH = Lambda.rows(obs_vars);

  arma::mat hh0;
  if (det(P0) == 0) {
   hh0 = h0;
  } else {
   hh0 = h0 + trans(chol(P0)) * as<arma::vec>(rnorm(iq));
  }

  mhh.row(0) = trans(mF * hh0 + mQ * trans(mu.row(0)));
  for (int i = 1; i < iT; i++) {
   mhh.row(i) = trans(mF * trans(mhh.row(i-1)) + mQ * trans(mu.row(i)));
  }

  /* The return is the smoothed state vector */
  return(mhh);
}

// [[Rcpp::export]]
arma::mat generate_mhh(       arma::mat mZ, arma::mat Lambda, arma::mat mF, arma::mat mQ, int iT, int ip, int iq, arma::mat h0, arma::mat P0) {
  /* This function generates the pseudo-state */
  /****************************************************/
  arma::mat mE(iT, iq);
  for (int i = 0; i < iq; i++) {
    mE.col(i) = as<arma::vec>(rnorm(iT));
  }

  arma::mat mZZ(iT, ip);
  mZZ.fill(NA_REAL);
  arma::mat mhh(iT, iq);
  mhh.fill(NA_REAL);

  arma::mat hh0;
  if (det(P0) == 0) {
    hh0 = h0;
  } else {
    hh0 = h0 + trans(chol(P0)) * as<arma::vec>(rnorm(iq));
  }
  mhh.row(0) = trans(mF * hh0 + mQ * trans(mE.row(0)));
  for (int i = 1; i < iT; i++) {
    mhh.row(i) = trans(mF * trans(mhh.row(i-1)) + mQ * trans(mE.row(i)));
  }

  return(mhh);
}

// [[Rcpp::export]]
arma::mat simulation_smoother(arma::mat mZ, arma::mat Lambda, arma::mat mF, arma::mat mQ, int iT, int ip, int iq, arma::mat h0, arma::mat P0) {
  /* This function produces a draw from the posterior distribution */
  /****************************************************/
  arma::mat mhh = generate_mhh(mZ, Lambda, mF, mQ, iT, ip, iq, h0, P0);
  arma::mat mZZ(iT, ip);
  arma::mat mz;
  arma::uvec obs_vars;
  arma::mat mH;
  arma::mat mzz;
  for (int i = 0; i < iT; i++) {
    mz = mZ.row(i);
    obs_vars = find_finite(mz);
    mH = Lambda.rows(obs_vars);
    mzz = mZZ.row(i);
    mzz.cols(obs_vars) = trans(mH * trans(mhh.row(i)));
    mZZ.row(i) = mzz;
  }

  arma::mat Z1 = smoother(mZ,  Lambda, mF, mQ, iT, ip, iq, h0, P0);
  arma::mat Z2 = smoother(mZZ, Lambda, mF, mQ, iT, ip, iq, h0, P0);
  return(Z1-Z2+mhh);
}
