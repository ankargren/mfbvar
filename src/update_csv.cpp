// Copyright of original code: Gregor Kastner (stochvol package)
// Copyright of modified code: Sebastian Ankargren (mfbvar package)
// The following code is a derivative work of the code
// developed by Gregor Kastner for the stochvol package, which
// is licensed GPL>=2. This code is therefore licensed under
// the terms of the GNU Public License, version 3.

#include "mfbvar.h"
#include "progutils.h"
#include "auxmix.h"
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
    const double prior_df) {
  // data: data matrix
  // phi: AR(1) parameter
  // sigma: standard deviation of log-volatility innovation
  // h: vector of log volatilities
  // h0: log volatility initial value
  // mixprob: mixture probabilities for Kim, Shephard, Chib (1998) algorithm
  // r: mixture indicators for KSC (1998) algorithm
  // priorlatent0: prior variance for initial value of log volatility
  // phi_invvar: inverse of prior variance for AR(1) parameter
  // phi_meaninvvar: prior mean of AR(1) parameter times phi_invvar
  // prior_sigma2: prior mean of variance of innovation
  // prior_df: prior degrees of freedom for variance of innovation
  int T = data.n_rows;
  int n = data.n_cols;

  arma::vec omega_diag(T+1);  // contains diagonal elements of precision matrix
  double omega_offdiag;  // contains off-diag element of precision matrix (const)
  arma::vec chol_offdiag(T), chol_diag(T+1);  // Cholesky-factor of Omega
  arma::vec covector(T+1);  // holds covector (see McCausland et al. 2011)
  arma::vec htmp(T+1);  // intermediate vector for sampling h
  arma::vec hnew(T+1);  // intermediate vector for sampling h

  double sigma2inv = std::pow(sigma, -2.0);

  double Bh0inv = 1.0/priorlatent0;

  arma::vec hT = h.rows(1, T - 1);
  arma::vec hT1 = h.rows(0, T - 2);


  /*
  * Sample phi
  */
  double phi_postvar = std::pow(phi_invvar + sigma2inv * arma::accu(arma::pow(hT1, 2.0)), -1.0);
  double phi_postmean = phi_postvar * (sigma2inv * arma::accu(hT1 % hT) + phi_meaninvvar);
  phi = rtruncnorm(phi_postmean, phi_postvar);
  const double phi2 = std::pow(phi, 2.0);

  /*
  * Sample sigma2
  */
  arma::vec u = hT - phi * hT1;
  sigma = std::pow(R::rgamma(prior_df + T - 1, 1/(prior_df * prior_sigma2 + arma::accu(arma::pow(u, 2.0)))), -0.5);
  sigma2inv = std::pow(sigma, -2.0);

  /*
  * Step (c): sample indicators
  */
  // calculate non-normalized CDF of posterior indicator probs
  for (int i = 0; i < n; ++i) {
    arma::vec mixprob_vec = arma::vec(10*T);
    arma::ivec r_vec = arma::ivec(T);
    findMixCDF(mixprob_vec, data.col(i)-h);
    invTransformSampling(mixprob_vec, r_vec, T);
    mixprob.col(i) = mixprob_vec;
    r.col(i) = r_vec;
  }


  // find correct indicators (currently by inversion method)

  /*
  * Step (a): sample the latent volatilities h:
  */
  omega_diag[0] = (Bh0inv + 1) * sigma2inv;
  covector[0] = 0.0;

  for (int j = 1; j < T; j++) {
    omega_diag[j] = sigma2inv*(1+phi2);
    covector[j] = 0.0;
    for (int i = 0; i < n; i++) {
      omega_diag[j] += mix_varinv[r.at(j-1, i)];
      covector[j] += (data.at(j-1, i) - mix_mean[r.at(j-1, i)])*mix_varinv[r.at(j-1, i)];
    }
  }
  omega_diag[T] = sigma2inv;
  covector[T] = 0.0;
  for (int i = 0; i < n; i++) {
    omega_diag[T] += mix_varinv[r.at(T-1, i)];
    covector[T] += (data.at(T-1, i) - mix_mean[r.at(T-1, i)])*mix_varinv[r.at(T-1, i)];
  }
  omega_offdiag = -phi * sigma2inv;  // omega_offdiag is constant

  // Cholesky decomposition
  cholTridiag(omega_diag, omega_offdiag, chol_diag, chol_offdiag);

  // Solution of Chol*x = covector ("forward algorithm")
  forwardAlg(chol_diag, chol_offdiag, covector, htmp);

  htmp += Rcpp::as<arma::vec>(Rcpp::rnorm(T+1));

  // Solution of (Chol')*x = htmp ("backward algorithm")
  backwardAlg(chol_diag, chol_offdiag, htmp, hnew);

  h = hnew.tail(T);
  h0 = hnew[0];
}
