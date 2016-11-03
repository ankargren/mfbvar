// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// build_U_cpp
mat build_U_cpp(mat Pi, int n_determ, int n_vars, int n_lags);
RcppExport SEXP MFBVAR_build_U_cpp(SEXP PiSEXP, SEXP n_determSEXP, SEXP n_varsSEXP, SEXP n_lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< int >::type n_determ(n_determSEXP);
    Rcpp::traits::input_parameter< int >::type n_vars(n_varsSEXP);
    Rcpp::traits::input_parameter< int >::type n_lags(n_lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_U_cpp(Pi, n_determ, n_vars, n_lags));
    return rcpp_result_gen;
END_RCPP
}
// build_demeaned_z_cpp
mat build_demeaned_z_cpp(mat z, mat psi, mat d, int n_T, int n_vars);
RcppExport SEXP MFBVAR_build_demeaned_z_cpp(SEXP zSEXP, SEXP psiSEXP, SEXP dSEXP, SEXP n_TSEXP, SEXP n_varsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< mat >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< mat >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n_T(n_TSEXP);
    Rcpp::traits::input_parameter< int >::type n_vars(n_varsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_demeaned_z_cpp(z, psi, d, n_T, n_vars));
    return rcpp_result_gen;
END_RCPP
}
// rmatnorm
mat rmatnorm(mat M, mat Q, mat P);
RcppExport SEXP MFBVAR_rmatnorm(SEXP MSEXP, SEXP QSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type M(MSEXP);
    Rcpp::traits::input_parameter< mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(rmatnorm(M, Q, P));
    return rcpp_result_gen;
END_RCPP
}
// rinvwish
mat rinvwish(int v, mat S);
RcppExport SEXP MFBVAR_rinvwish(SEXP vSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvwish(v, S));
    return rcpp_result_gen;
END_RCPP
}
// rmultn
vec rmultn(vec m, mat Sigma);
RcppExport SEXP MFBVAR_rmultn(SEXP mSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultn(m, Sigma));
    return rcpp_result_gen;
END_RCPP
}
