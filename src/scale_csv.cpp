#include "mfbvar.h"
void scale_csv(arma::vec & exp_sqrt_f, arma::mat & y_scaled, arma::mat & X_scaled,
          const arma::mat & y_in, const arma::mat & X_in,
          const arma::vec & f_in) {
  exp_sqrt_f = arma::exp(0.5 * f_in);
  y_scaled = y_in;
  y_scaled.each_col() /= exp_sqrt_f;
  X_scaled = X_in;
  X_scaled.each_col() /= exp_sqrt_f;
}
