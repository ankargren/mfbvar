#ifndef MFBVAR_MFBVAR_H
#define MFBVAR_MFBVAR_H

#include <RcppArmadillo.h>
#include "mvn.h"
#include "simsm_adaptive_univariate.h"
#include "simsm_adaptive_cv.h"

arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags);
arma::vec rmultn(arma::vec m, arma::mat Sigma);

#endif
