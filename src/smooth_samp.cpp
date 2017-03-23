// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h> 

using namespace arma; 
using namespace Rcpp;

vec vec_na_rm(const vec vx){
  vec vr(vx.size());
  int k = 0;
  for(int i=0; i<vx.size(); i++){
    if(is_finite(vx(i))) {
      vr(k)=vx(i); k++;
    }
  }
  vr.resize(k);
  return(vr);
}



// [[Rcpp::export]]

mat smooth(const mat & mZ, const mat & mX, const List & lH, mat & mF, const mat & mB, const mat & mQ,const int iT, const int ip, const int iq, const int is, const vec & h0, const mat & P0, const vec & X0){

  // initialize 
  mat QQ = mQ * mQ.t(); // q by q
  List mv(iT), IS(iT), aK(iT);

  //// filtering 

  // Predict
  vec h1 = mF * h0 + mB * X0; // q by 1
  mat P1 = mF * P0 * mF.t() + QQ; // q by q
  // Update
  mat mH = as<mat>(lH[0]); // p by q
  vec vz = vec_na_rm(vectorise(mZ.row(0)));
  
  mat mvv = vz - mH * h1; // smoothing
  mv[0] = mvv;
  mat aS = mH * P1 * mH.t();
  mat iS = inv_sympd(aS); // smoothing
  IS[0] = iS; 
  mat akk = P1 * mH.t() * iS; // smoothing
  aK[0] = akk;
  vec h2 = h1 + akk * mvv;
  mat idq = eye<mat>(iq,iq);
  mat P2 = (idq - akk * mH) * P1;

  int iter;
  for(iter=1; iter<iT; iter++){
    // Predict
    h1 = mF * h2 + mB * vectorise(mX.row(iter-1));
    P1 = mF * P2 * mF.t() + QQ;
    // Update
    mH = as<mat>(lH[iter]);
    vz = vec_na_rm(vectorise(mZ.row(iter)));
    
    mvv = vz - mH * h1; // smoothing 
    mv[iter] = mvv;
    aS = mH * P1 * mH.t();
    iS = inv_sympd(aS); // smoothing
    IS[iter] = iS;
    akk = P1 * mH.t() * iS; // smoothing
    aK[iter] = akk;
    h2 = h1 + akk * mvv;
    P2 = (idq - akk * mH) * P1; 
  }

  // smoothing
  
  mat mr(iT,iq,fill::zeros), mu(iT,iq), fk, mee;

  for(iter=iT-1; iter>0; iter--){
    mH = as<mat>(lH[iter]);
    akk = as<mat>(aK[iter]);
    fk = mF * akk;

    iS = as<mat>(IS[iter]);
    mvv = as<mat>(mv[iter]);

    mee = iS * mvv - (mr.row(iter) * fk).t();
    mr.row(iter-1) = (mH.t() * mee + mF.t() * mr.row(iter).t()).t();
    mu.row(iter) = mr.row(iter) * mQ;
  }

  mH = as<mat>(lH[0]);
  akk = as<mat>(aK[0]);
  fk = mF * akk;

  iS = as<mat>(IS[0]);
  mvv = as<mat>(mv[0]);

  mee = iS * mvv - (mr.row(0)*fk).t();
  mat r0 = mH.t() * mee + mF.t() * mr.row(0).t();
  mu.row(0) = mr.row(0) * mQ;

  mu = join_cols((mQ.t() * r0).t(), mu);

  return(mu); 
}   
