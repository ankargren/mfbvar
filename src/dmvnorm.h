#ifndef _DMVNORM_H
#define _DMVNORM_H

//#define ARMA_NO_DEBUG // disables bounds checks
#include <RcppArmadillo.h>

// Main predict function (as called from R):
RcppExport SEXP dmvnorm(const SEXP, const SEXP, const SEXP, const SEXP);

#endif
