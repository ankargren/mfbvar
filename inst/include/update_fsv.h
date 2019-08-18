#ifndef _UPDATE_FSV_H_
#define _UPDATE_FSV_H_
void update_fsv(arma::mat & armafacload, arma::mat & armaf, arma::mat & armah,
                arma::vec & armah0,
                Rcpp::NumericMatrix & curpara,
                const arma::mat & armatau2,
                const arma::mat & armay,
                const double bmu, const double Bmu, const double a0idi, const double b0idi,
                const double a0fac, const double b0fac, const Rcpp::NumericVector & Bsigma,
                const double B011inv, const double B022inv,
                const Rcpp::NumericVector & priorh0, const arma::imat & armarestr);
#endif
