#include "mfbvar.h"

/*
 double gigrvg(double lambda, double chi, double psi) {
 SEXP gig = do_rgig(1, lambda, chi, psi);
 double ret = 0.0;//Rcpp::as<double>(gig);
 return ret;
 }
 */

// [[Rcpp::export]]
double do_rgig1(double lambda, double chi, double psi) {
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  return Rcpp::as<double>(fun(1, lambda, chi, psi));
}

double rig(double mu, double lambda){
  double z = R::rnorm(0,1);
  double y = z*z;
  double x = mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
  double u=R::runif(0,1);
  double out;
  if(u <= mu/(mu+x)){
    out = x;
  } else {
    out = mu*mu/x;
  }
  return out;
}
