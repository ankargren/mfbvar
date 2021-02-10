#include "progutils_fsv.h"
// Copyright of original code: Gregor Kastner (factorstochvol package)
// Copyright of modified code: Sebastian Ankargren (mfbvar package)
// The following code is a derivative work of the code
// developed by Gregor Kastner for the factorstochvol package, which
// is licensed GPL>=2. This code is therefore licensed under
// the terms of the GNU Public License, version 3.

double logdnormquot(double x, double y, double mu, double sigma) {
 return ((y-mu)*(y-mu) - (x-mu)*(x-mu)) / (2*sigma*sigma);
}

double logspecialquot(double x, double y, double alpha, double beta, double c) {
 return (alpha/c) * (x - y) - beta * (exp(x/c) - exp(y/c));
}


void store(const Rcpp::NumericMatrix &curfacload, Rcpp::NumericVector &facload,
           const Rcpp::NumericMatrix &curf,       Rcpp::NumericVector &f,
	   const Rcpp::NumericMatrix &curh,       Rcpp::NumericVector &h,
	   const Rcpp::NumericVector &curh0,      Rcpp::NumericMatrix &h0,
	   const Rcpp::NumericMatrix &curpara,    Rcpp::NumericVector &para,
	   const Rcpp::NumericVector &curlambda2, Rcpp::NumericMatrix &lambda2,
	   const Rcpp::NumericMatrix &curtau2,    Rcpp::NumericVector &tau2,
           const Rcpp::NumericVector &curmixprob, Rcpp::NumericVector &mixprob,
	   const Rcpp::IntegerMatrix &curmixind,  Rcpp::IntegerVector &mixind,
	   const bool auxstore, const int thintime, const int where) {

 std::copy(curfacload.begin(), curfacload.end(), facload.begin() + where * curfacload.length());
 std::copy(curpara.begin(), curpara.end(), para.begin() + where * curpara.length());

 if (thintime == 1) { // store everything

  std::copy(curf.begin(), curf.end(), f.begin() + where * curf.length());
  std::copy(curh.begin(), curh.end(), h.begin() + where * curh.length());

 } else if (thintime == -1) { // store only t = T

  for (int i = 0; i < curf.nrow(); i++) {
   f(where*curf.nrow() + i) = curf(i, curf.ncol()-1);
  }

  for (int i = 0; i < curh.ncol(); i++) {
   h(where*curh.ncol() + i) = curh(curh.nrow()-1, i);
  }

 } else if (thintime > 1) { // store every thintimeth point in time

  int tmp = curf.ncol()/thintime;
  int tmpp = where * curf.nrow() * tmp;

  for (int j = 0; j < tmp; ++j) {
   int tmppp = j*thintime;
   int tmpppp = tmpp + j*curf.nrow();

   for (int i = 0; i < curf.nrow(); ++i) {
    f(tmpppp + i) = curf(i, tmppp);
   }
  }

  tmpp = where * curh.ncol() * tmp;

  for (int i = 0; i < curh.ncol(); ++i) {
   int tmpppp = tmpp + i*tmp;

   for (int j = 0; j < tmp; ++j) {
    h(tmpppp + j) = curh(j*thintime, i);
   }
  }
 }

 std::copy(curh0.begin(), curh0.end(), h0.begin() + where * curh0.length());

 if (auxstore) { // store mixture probabilities, mixture indicators, shrinkage hyperparas, h0
  std::copy(curmixprob.begin(), curmixprob.end(), mixprob.begin() + where * curmixprob.length());
  std::copy(curmixind.begin(), curmixind.end(), mixind.begin() + where * curmixind.length());
  std::copy(curlambda2.begin(), curlambda2.end(), lambda2.begin() + where * curlambda2.length());
  std::copy(curtau2.begin(), curtau2.end(), tau2.begin() + where * curtau2.length());
 }
}
