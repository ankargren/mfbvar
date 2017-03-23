// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#define _USE_MATH_DEFINES // for C++
#include <cmath>
// [[Rcpp::export]]
arma::vec vec_na_rm(const arma::vec vx){
  arma::vec vr(vx.size());
  uword k = 0;
  for(uword i=0; i<vx.size(); i++){
    if(is_finite(vx(i))) {
      vr(k)=vx(i); k++;
    }
  }
  vr.resize(k);
  return(vr);
}



// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#define _USE_MATH_DEFINES // for C++
#include <cmath>
// [[Rcpp::export]]
arma::mat smooth(const arma::mat & mZ, const arma::mat & mX, const Rcpp::List & lH, arma::mat & mF, const arma::mat & mB, const arma::mat & mQ,const int iT, const int ip, const int iq, const int is, const arma::vec & h0, const arma::mat & P0, const arma::vec & X0){

  // initialize
  arma::mat QQ = mQ * mQ.t(); // q by q
  Rcpp::List mv(iT), IS(iT), aK(iT);

  //// filtering

  // Predict
  arma::vec h1 = mF * h0 + mB * X0; // q by 1
  arma::mat P1 = mF * P0 * mF.t() + QQ; // q by q
  // Update
  arma::mat mH = as<arma::mat>(lH[0]); // p by q
  arma::vec vz = vec_na_rm(vectorise(mZ.row(0)));

  arma::mat mvv = vz - mH * h1; // smoothing
  mv[0] = mvv;
  arma::mat aS = mH * P1 * mH.t();
  arma::mat iS = arma::inv_sympd(aS); // smoothing
  IS[0] = iS;
  arma::mat akk = P1 * mH.t() * iS; // smoothing
  aK[0] = akk;
  arma::vec h2 = h1 + akk * mvv;
  arma::mat idq = eye<arma::mat>(iq,iq);
  arma::mat P2 = (idq - akk * mH) * P1;

  int iter;
  for(iter=1; iter<iT; iter++){
    // Predict
    h1 = mF * h2 + mB * vectorise(mX.row(iter-1));
    P1 = mF * P2 * mF.t() + QQ;
    // Update
    mH = as<arma::mat>(lH[iter]);
    vz = vec_na_rm(vectorise(mZ.row(iter)));

    mvv = vz - mH * h1; // smoothing
    mv[iter] = mvv;
    aS = mH * P1 * mH.t();
    iS = arma::inv_sympd(aS); // smoothing
    IS[iter] = iS;
    akk = P1 * mH.t() * iS; // smoothing
    aK[iter] = akk;
    h2 = h1 + akk * mvv;
    P2 = (idq - akk * mH) * P1;
  }

  // smoothing

  arma::mat mr(iT,iq,fill::zeros), mu(iT,iq), fk, mee;

  for(iter=iT-1; iter>0; iter--){
    mH = as<arma::mat>(lH[iter]);
    akk = as<arma::mat>(aK[iter]);
    fk = mF * akk;

    iS = as<arma::mat>(IS[iter]);
    mvv = as<arma::mat>(mv[iter]);

    mee = iS * mvv - (mr.row(iter) * fk).t();
    mr.row(iter-1) = (mH.t() * mee + mF.t() * mr.row(iter).t()).t();
    mu.row(iter) = mr.row(iter) * mQ;
  }

  mH = as<arma::mat>(lH[0]);
  akk = as<arma::mat>(aK[0]);
  fk = mF * akk;

  iS = as<arma::mat>(IS[0]);
  mvv = as<arma::mat>(mv[0]);

  mee = iS * mvv - (mr.row(0)*fk).t();
  arma::mat r0 = mH.t() * mee + mF.t() * mr.row(0).t();
  mu.row(0) = mr.row(0) * mQ;

  mu = arma::join_cols((mQ.t() * r0).t(), mu);

  return(mu);
}
