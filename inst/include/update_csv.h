#ifndef _UPDATE_CSV_H_
#define _UPDATE_CSV_H_
void update_csv(
    const arma::mat& data,
    double& phi,
    double& sigma,
    arma::vec& h,
    double& h0,
    arma::mat& mixprob,
    arma::imat& r,
    const double priorlatent0,
    const double phi_invvar,
    const double phi_meaninvvar,
    const double prior_sigma2,
    const double prior_df,
    bool fixate_latent,
    bool fixate_latent0,
    bool fixate_phi,
    bool fixate_sigma);

void scale_csv(arma::vec & exp_sqrt_f, arma::mat & y_scaled,
               arma::mat & X_scaled, const arma::mat & y_in,
               const arma::mat & X_in, const arma::vec & f_in);

void fcst_csv(arma::mat & Z_fcst, const arma::mat & Z, const arma::mat & Pi,
              const arma::mat & Sigma_chol, double phi, double sigma, double vol_pred,
              arma::uword n_fcst, arma::uword n_lags, arma::uword n_vars);

#endif
