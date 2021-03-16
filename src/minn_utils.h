#ifndef _MINN_UTILS_H_
#define _MINN_UTILS_H_
arma::mat create_X_t(const arma::mat & y);
arma::mat create_X(const arma::mat & y, arma::uword k);
arma::mat sample_Pi_mat(const arma::mat & post_Pi,
                   const arma::mat & post_Pi_Omega,
                   const arma::mat & Sigma,
                   bool check_roots,
                   arma::uword n_vars);
#endif
