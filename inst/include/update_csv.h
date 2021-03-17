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

#endif

