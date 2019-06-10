#include <RcppParallel.h>
#include "mvn.h"
struct Pi_parallel : public RcppParallel::Worker {
  arma::mat & output;
  const arma::mat & y;
  const arma::mat & X;
  const arma::mat & d;
  const arma::mat & eps;
  const arma::mat & factors;
  arma::uword T;
  arma::uword n;
  arma::uword p;

  Pi_parallel(arma::mat & output,
              const arma::mat & y,
              const arma::mat & X,
              const arma::mat & d,
              const arma::mat & eps,
              const arma::mat & factors,
              const arma::uword T,
              const arma::uword n,
              const arma::uword p);

  void operator()(std::size_t begin, std::size_t end);
};

arma::mat sample_Pi_parallel(const arma::mat & y, const arma::mat & d, const arma::mat & eps,
                             const arma::mat & factors,
                             const int T, const int n, const int p);
