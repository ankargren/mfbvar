#include <RcppParallel.h>
#include "mvn.h"
struct Pi_parallel_rue : public RcppParallel::Worker {
  arma::mat & output;
  const arma::mat & y;
  const arma::mat & X;
  const arma::mat & d;
  const arma::mat & eps;
  const arma::mat & volatility;
  const arma::mat & prior_AR1;
  arma::uword T;
  arma::uword n;
  arma::uword p;

  Pi_parallel_rue(arma::mat & output,
              const arma::mat & y,
              const arma::mat & X,
              const arma::mat & d,
              const arma::mat & eps,
              const arma::mat & volatility,
              const arma::mat & prior_AR1,
              const arma::uword T,
              const arma::uword n,
              const arma::uword p);

  void operator()(std::size_t begin, std::size_t end);
};

struct Pi_parallel_bcm : public RcppParallel::Worker {
  arma::mat & output;
  const arma::mat & y;
  const arma::mat & X;
  const arma::mat & d;
  const arma::mat & eps;
  const arma::mat & volatility;
  arma::uword T;
  arma::uword n;
  arma::uword p;

  Pi_parallel_bcm(arma::mat & output,
                  const arma::mat & y,
                  const arma::mat & X,
                  const arma::mat & d,
                  const arma::mat & eps,
                  const arma::mat & volatility,
                  const arma::uword T,
                  const arma::uword n,
                  const arma::uword p);

  void operator()(std::size_t begin, std::size_t end);
};
