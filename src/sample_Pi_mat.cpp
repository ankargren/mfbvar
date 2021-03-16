#include "mfbvar.h"
#include "ss_utils.h"

arma::mat sample_Pi_mat(const arma::mat & post_Pi,
                        const arma::mat & post_Pi_Omega,
                        const arma::mat & Sigma,
                        bool check_roots,
                        arma::uword n_vars) {
  bool stationarity_check = false;
  int num_try = 0, iter = 0;
  double root = 1000;
  arma::mat Pi, Pi_comp;
  while (stationarity_check == false) {
    iter += 1;
    Pi = rmatn(post_Pi.t(), post_Pi_Omega, Sigma);
    if (check_roots) {
      Pi_comp.rows(0, n_vars-1) = Pi;
      root = max_eig_cpp(Pi_comp);
    } else {
      root = 0.0;
    }
    if (root < 1.0) {
      stationarity_check = true;
      num_try = iter;
    }
    if (iter == 1000) {
      Rcpp::stop("Attemped to draw stationary Pi 1,000 times.");
    }
  }
  return Pi;
}
