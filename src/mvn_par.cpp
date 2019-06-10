#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "mvn.h"
#include "mvn_par.h"

Pi_parallel::Pi_parallel(arma::mat & output,
              const arma::mat & y,
              const arma::mat & X,
              const arma::mat & d,
              const arma::mat & eps,
              const arma::mat & factors,
              const arma::uword T,
              const arma::uword n,
              const arma::uword p) : output(output), y(y), X(X), d(d), eps(eps), factors(factors), T(T), n(n), p(p) {};

void Pi_parallel::operator()(std::size_t begin, std::size_t end) {
  for (std::size_t i = begin; i < end; i++)
  {
    arma::vec y_i = arma::vec(y.begin_col(i)+p, T);
    arma::vec factors_i = arma::exp(-0.5*factors.col(i));
    y_i %= factors_i;
    arma::mat X_i = X.each_col() % factors_i;
    arma::vec eps_i = eps.unsafe_col(i);
    arma::vec d_i = d.unsafe_col(i);
    arma::vec armatmp = output.unsafe_col(i);
    armatmp = mvn_rue_eps(X_i, d_i, y_i, eps_i);
  }
}

// [[Rcpp::export]]
arma::mat sample_Pi_parallel(const arma::mat & y, const arma::mat & d, const arma::mat & eps,
                             const arma::mat & factors,
                             const int T, const int n, const int p) {
  arma::mat output(n*p+1, n);
  arma::mat X = arma::mat(T-p, n*p+1);
  X.col(0).fill(1.0);
  for (arma::uword j = 0; j<p; ++j) {
    X.cols(1+0+j*n, 1+(j+1)*n-1) = y.rows(p-j-1, T-2-j);
  }
  Pi_parallel Pi_parallel_i(output, y, X, d, eps, factors, T, n, p);
  RcppParallel::parallelFor(0, n, Pi_parallel_i);
  return output.t();
}

