// Functions related to sampling the indicators
// Constants relating to the approximation of log(chisq) through
// normal mixture (Omori et al., 2007) can be found in auxmix.h

#include <RcppArmadillo.h>
#include "auxmix.h"

using namespace Rcpp;

// Non-normalized posterior probabilities
void findMixprobs(
    arma::vec& mixprob,
    const arma::vec& datanorm)  {
 int T = datanorm.size();
 int tmp;
 for (int c = 0; c < T; c++) {  // SLOW (10*T calls to exp)!
  tmp = 10*c;
  for (int r = 0; r < 10; r++) {
   mixprob[tmp+r] = exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
  }
 }
}

// Cumulative sum over columns of a matrix
void colCumsums(
    arma::vec& x,
    int const nrow,
    int const ncol) {
 int tmp;
 for (int c = 0; c < ncol; c++) {
  tmp = c*nrow;
  for (int r = 1; r < nrow; r++) {
   x[tmp+r] = x[tmp+r-1] + x[tmp+r];
  }
 }
}

// Combines findMixprobs() and colCumsums() (see above) into one function
void findMixCDF(
    arma::vec& mixprob,
    const arma::vec& datanorm)  {
 int T = datanorm.size();
 int tmp;
 for (int c = 0; c < T; c++) {  // SLOW (10*T calls to exp)!
  tmp = 10*c;
  mixprob[tmp] = exp(mix_pre[0]-(datanorm[c]-mix_mean[0])*(datanorm[c]-mix_mean[0])*mix_2varinv[0]);
  for (int r = 1; r < 10; r++) {
   mixprob[tmp+r] = mixprob[tmp+r-1] + exp(mix_pre[r]-(datanorm[c]-mix_mean[r])*(datanorm[c]-mix_mean[r])*mix_2varinv[r]);
  }
 }
}

void invTransformSampling(
    const arma::vec& mixprob,
    arma::ivec& r,
    int T) {
  int index;
  arma::vec innov = runif(T);
  double temp;
  bool larger, smaller;
  for (int j = 0; j < T; j++) {
    index = (10-1)/2;  // start searching in the middle
    temp = innov[j]*mixprob[9 + 10*j];  // current (non-normalized) value
    larger = false;  // indicates that we already went up
    smaller = false; // indicates that we already went down
    while(true) {
      if (temp > mixprob[index +  10*j]) {
        if (smaller == true) {
          index++;
          break;
        }
        else {
          index++;
          larger = true;
        }
      }
      else {
        if (larger == true) {
          break;
        }
        else {
          if (index == 0) {
            break;
          }
          else {
            index--;
            smaller = true;
          }
        }
      }
    }
    r[j] = index;
  }
}
