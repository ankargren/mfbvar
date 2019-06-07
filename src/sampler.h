#ifndef _SAMPLER_H
#define _SAMPLER_H

//#define ARMA_NO_DEBUG // disables bounds checks
#include <RcppArmadillo.h>
#include <stochvol.h>  // decl'd and def'd in "stochvol" (univariate SV-update)
#include "progutils.h"

double rgig1(double, double, double);

// Main sampler (as called from R):
RcppExport SEXP sampler(const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
  const SEXP, const SEXP, const SEXP, const SEXP);

RcppExport SEXP sampler2(const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP, const SEXP,
                        const SEXP, const SEXP, const SEXP, const SEXP);

#endif
