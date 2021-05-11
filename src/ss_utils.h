#ifndef _SS_UTILS_H_
#define _SS_UTILS_H_
arma::mat build_U_cpp(const arma::mat & Pi, int n_determ, int n_vars, int n_lags);
void update_demean(arma::mat & my, arma::mat & mu_long,
                   const arma::mat & y_in_p, const arma::mat & mu_mat, const arma::mat & d1,
                   const arma::mat & Psi_i, const arma::mat & Lambda_single,
                   arma::uword n_vars, arma::uword n_q, arma::uword n_Lambda, arma::uword n_T);

void posterior_psi_iw(arma::vec & psi_i, arma::mat & mu_mat,
                      const arma::mat & Pi_i, const arma::mat & D_mat,
                      const arma::mat & Sigma_i, const arma::mat & inv_prior_psi_Omega,
                      const arma::mat & Z_i, const arma::mat & X,
                      const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                      int n_determ, int n_vars, int n_lags);
void posterior_psi_csv(arma::vec & psi_i, arma::mat & mu_mat,
                       const arma::mat & Pi_i, const arma::mat & D_mat,
                       const arma::mat & Sigma_chol_inv, const arma::mat & exp_sqrt_f,
                       const arma::mat & inv_prior_psi_Omega,
                       const arma::mat & Z_i, const arma::mat & X,
                       const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                       int n_determ, int n_vars, int n_lags);
void posterior_psi_fsv(arma::vec & psi_i, arma::mat & mu_mat,
                       const arma::mat & Pi_i, const arma::mat & D_mat,
                       const arma::mat & idivar, const arma::mat & inv_prior_psi_Omega,
                       const arma::mat & Z_i, const arma::mat & X,
                       const arma::mat & startfacload, const arma::mat & startfac,
                       const arma::mat & inv_prior_psi_Omega_mean, const arma::mat & dt,
                       int n_determ, int n_vars, int n_lags);

double max_eig_cpp(const arma::mat & A);
#endif
