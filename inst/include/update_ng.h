#ifndef _UPDATE_NG_H_
#define _UPDATE_NG_H_
void update_ng(double & phi_mu, double & lambda_mu, arma::vec & omega,
               arma::uword nm, const double c0, const double c1, double s,
               const arma::vec & psi_i, const arma::vec & prior_psi_mean,
               arma::uword & accept, bool fixate_phi_mu, bool fixate_lambda_mu,
               bool fixate_omega);
void update_s(double & s,
              arma::running_stat<double> & stats,
              double accept,
              arma::uword i,
              double M);
#endif
