#ifndef _PROGUTILS_H
#define _PROGUTILS_H

#include <RcppArmadillo.h>

double logdnormquot(double x, double y, double mu, double sigma);
double logspecialquot(double x, double y, double alpha, double beta, double c);

void store(const Rcpp::NumericMatrix &curfacload, Rcpp::NumericVector &facload,
           const Rcpp::NumericMatrix &curf,       Rcpp::NumericVector &f,
	   const Rcpp::NumericMatrix &curh,       Rcpp::NumericVector &h,
	   const Rcpp::NumericVector &curh0,      Rcpp::NumericMatrix &h0,
	   const Rcpp::NumericMatrix &curpara,    Rcpp::NumericVector &para,
	   const Rcpp::NumericVector &curlambda2, Rcpp::NumericMatrix &lambda2,
	   const Rcpp::NumericMatrix &curtau2,    Rcpp::NumericVector &tau2,
           const Rcpp::NumericVector &curmixprob, Rcpp::NumericVector &mixprob,
	   const Rcpp::IntegerMatrix &curmixind,  Rcpp::IntegerVector &mixind,
	   const bool auxstore, const int thintime, const int where);

#endif
