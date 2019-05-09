#ifndef _SS_UTILS_H_
#define _SS_UTILS_H_
arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags);
arma::mat create_X_t_noint(const arma::mat & y);
arma::mat create_X_noint(const arma::mat & y, arma::uword k);
#endif
