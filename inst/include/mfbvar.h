#ifndef MFBVAR_MFBVAR_H
#define MFBVAR_MFBVAR_H

#include <RcppArmadillo.h>
#include "mvn.h"
#include "simsm_adaptive_univariate.h"
#include "simsm_adaptive_cv.h"
#include "simsm_adaptive_sv.h"
#include <progress.hpp>
#include "progress_bar.hpp"



arma::vec rmultn(const arma::vec & m, const arma::mat & Sigma);
arma::mat rinvwish(int v, const arma::mat & S);
arma::mat rmatn(const arma::mat & M, const arma::mat & Q, const arma::mat & P);

// Import the rgig
double do_rgig1(double lambda, double chi, double psi);
double rig(double mu, double lambda);

#endif
