#ifndef _SS_UTILS_H_
#define _SS_UTILS_H_
arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags);
arma::mat create_X_t_noint(const arma::mat & y);
arma::mat create_X_noint(const arma::mat & y, arma::uword k);
void update_demean(arma::mat & my, arma::mat & mu_long,
                   const arma::mat & y_in_p, const arma::mat & mu_mat, const arma::mat & d1,
                   const arma::mat & Psi_i, const arma::mat & Lambda_single,
                   arma::uword n_vars, arma::uword n_q, arma::uword n_Lambda, arma::uword n_T);
#endif
