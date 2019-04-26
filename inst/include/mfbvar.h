#ifndef MFBVAR_MFBVAR_H
#define MFBVAR_MFBVAR_H

#include <RcppArmadillo.h>
#include "mvn.h"
#include "simsm_adaptive_univariate.h"
#include "simsm_adaptive_cv.h"

arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags);
arma::mat create_X_t(const arma::mat & y);
arma::mat create_X(const arma::mat & y, arma::uword k);
arma::mat create_X_t_noint(const arma::mat & y);
arma::mat create_X_noint(const arma::mat & y, arma::uword k);

arma::vec rmultn(const arma::vec & m, const arma::mat & Sigma);
arma::mat rinvwish(int v, const arma::mat & S);
arma::mat rmatn(const arma::mat & M, const arma::mat & Q, const arma::mat & P);

void posterior_psi_iw(arma::vec & psi_i, arma::mat & mu_mat,
                 const arma::mat & Pi_i, const arma::mat & D_mat,
                 const arma::mat & Sigma_i, const arma::mat & inv_prior_psi_Omega,
                 const arma::mat & Z_i, const arma::mat & X,
                 const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                 int n_determ, int n_vars, int n_lags);

double max_eig_cpp(const arma::mat & A);
#endif
