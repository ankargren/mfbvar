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
