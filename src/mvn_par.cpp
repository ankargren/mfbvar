#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "mvn.h"
#include "mvn_par.h"

Pi_parallel_rue::Pi_parallel_rue(arma::mat & output,
              const arma::mat & y,
              const arma::mat & X,
              const arma::mat & d,
              const arma::mat & eps,
              const arma::mat & volatility,
              const arma::mat & prior_AR1,
              const arma::uword T,
              const arma::uword n,
              const arma::uword p) : output(output), y(y), X(X), d(d), eps(eps), volatility(volatility), prior_AR1(prior_AR1), T(T), n(n), p(p) {};

void Pi_parallel_rue::operator()(std::size_t begin, std::size_t end) {
  for (std::size_t i = begin; i < end; i++)
  {
    arma::vec h_j = arma::exp(-0.5 * volatility.col(i));
    arma::mat X_j = X.each_col() % h_j;
    arma::vec y_j = y.col(i) % h_j;
    arma::vec eps_i = eps.unsafe_col(i);
    arma::vec d_i = d.unsafe_col(i);
    output.col(i) = mvn_rue_eps(X_j, d_i, y_j, eps_i, prior_AR1(i), i);
  }
}

Pi_parallel_bcm::Pi_parallel_bcm(arma::mat & output,
                                 const arma::mat & y,
                                 const arma::mat & X,
                                 const arma::mat & d,
                                 const arma::mat & eps,
                                 const arma::mat & volatility,
                                 const arma::uword T,
                                 const arma::uword n,
                                 const arma::uword p) : output(output), y(y), X(X), d(d), eps(eps), volatility(volatility), T(T), n(n), p(p) {};

void Pi_parallel_bcm::operator()(std::size_t begin, std::size_t end) {
  for (std::size_t i = begin; i < end; i++)
  {
    arma::vec h_j = arma::exp(-0.5 * volatility.col(i));
    arma::mat X_j = X.each_col() % h_j;
    arma::vec y_j = y.col(i) % h_j;
    arma::vec eps_i = eps.unsafe_col(i);
    arma::vec d_i = d.unsafe_col(i);
    output.col(i) = mvn_bcm_eps(X_j, d_i, y_j, eps_i);
  }
}


