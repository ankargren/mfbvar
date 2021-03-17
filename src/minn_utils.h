#ifndef _MINN_UTILS_H_
#define _MINN_UTILS_H_
arma::mat create_X_t(const arma::mat & y);
arma::mat create_X(const arma::mat & y, arma::uword k);
arma::mat sample_Pi_mat(const arma::mat & post_Pi,
                   const arma::mat & post_Pi_Omega,
                   const arma::mat & Sigma,
                   bool check_roots,
                   arma::uword n_vars);
arma::mat sample_Pi_vec(const arma::mat & Sigma_inv,
              const arma::mat & X,
              const arma::mat & y,
              const arma::vec & prior_Pi_Omega_vec_inv,
              const arma::mat & Omega_Pi,
              bool check_roots,
              arma::uword n_vars,
              arma::uword n_lags);
#endif
