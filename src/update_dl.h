#ifndef MFBVAR_UPDATE_DL_H
#define MFBVAR_UPDATE_DL_H
void update_dl(arma::mat & prior_Pi_Omega, arma::vec & aux,
               arma::vec & local, double & global, const arma::mat & Pi_i,
               arma::uword n_vars, arma::uword n_lags, const double a);
#endif
