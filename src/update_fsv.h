#include <RcppArmadillo.h>
#include <stochvol.h>  // decl'd and def'd in "stochvol" (univariate SV-update)
#include "progutils_fsv.h"
#include "auxmix.h"


// curfacload changed to armafacload
void update_fsv(arma::mat & armafacload, arma::mat & armaf, arma::mat & armah,
                arma::vec & armah0,
                Rcpp::NumericMatrix & curpara,
                const arma::mat & armatau2,
                const arma::mat & armay,
                const double bmu, const double Bmu, const double a0idi, const double b0idi,
                const double a0fac, const double b0fac, const Rcpp::NumericVector & Bsigma,
                  const double B011inv, const double B022inv, const bool Gammaprior,
                  const bool truncnormal, const double MHcontrol, const int MHsteps,
                  const int parameterization, const Rcpp::NumericVector & sv,
                  const Rcpp::NumericVector & priorhomoskedastic,
                  const Rcpp::NumericVector & priorh0, const arma::imat & armarestr) {



  const int interweaving       = 4;
  const bool signswitch        = true;
  const bool samplefac        = true;
  int nlambda = 0;
  const double c0 = 2.5;
  const Rcpp::NumericVector C0 = 1.5*Bsigma;

  using namespace Rcpp;

  const int m = armay.n_rows; // number of time series
  const int T = armay.n_cols; // length of time series
  const int r = armafacload.n_cols; // number of latent factors
  const int mpr = m + r;

  arma::irowvec nonzerospercol = arma::sum(armarestr, 0);
  arma::icolvec nonzerosperrow = arma::sum(armarestr, 1);

  // restriction on factor loadings matrix:
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < r; j++) {
      if (armarestr(i, j) == 0) armafacload(i,j) = 0.;
    }
  }
  // restriction on factor loadings matrix:
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < r; j++) {
      if (armarestr(i,j) == 0) armatau2(i,j) = 0.;
    }
  }
  // pre-calculation of a posterior parameter
  double cT = 0;
  if (Gammaprior) {
    if (MHsteps == 2 || MHsteps == 3) cT = T/2.0; // we want IG(-.5,0) as proposal
    else if (MHsteps == 1) cT = (T-1)/2.0; // we want IG(-.5,0) as proposal
  } else {
    if (MHsteps == 2) cT = c0 + (T+1)/2.0;  // pre-calculation outside the loop
  }

  int tmpcounter = 0;
  arma::uvec diagindices(m);
  for (int k = 0; k < m; k++) {
    for (int l = k; l < m; l++) {
      if (k == l) diagindices(k) = tmpcounter;
      tmpcounter++;
    }
  }
  for (int j = m; j < mpr; j++) {
    if (sv(j) == false) {
      armah.col(j).fill(0.);
    }
  }

  //convention: "arma"-prefixed variables denote Armadillo proxy objects
  arma::mat armafacloadt = arma::trans(armafacload);
  arma::uvec armafacloadtunrestrictedelements = arma::find(armarestr.t() != 0);
  arma::vec armafacloadtmp = arma::zeros<arma::vec>(armafacloadtunrestrictedelements.size());
  arma::vec armafacload2inter(r, arma::fill::zeros);
  arma::mat armahtilde(armah.n_rows, m);





  //current shrinkage latents lambda^2
  arma::vec armalambda2(nlambda);

  // temporary stroage for hopen in interweaving
  arma::vec hopen(T);

  // NOTE: (Almost) all storage of MCMC draws is done in NumericVectors
  // because no 'array' structure is available at this point in time.

  // facload holds the factor loadings:
  NumericVector facload(m * r);
  facload.attr("dim") = Dimension(m, r, 1);

  //current mixture indicator draws
  IntegerMatrix curmixind(T, mpr);
  NumericVector curmixprob(10 * T * mpr);

  // h holds the latent log-volatilities, but not h0!
  NumericVector h(T);
  h.attr("dim") = Dimension(1, T, 1);

  // f holds the latent factor draws
  NumericVector f(T);
  f.attr("dim") = Dimension(T, 1, 1);


  // mixind holds the mixture indicators for the auxiliary mixture sampling
  IntegerVector mixind(T * mpr);

  // mixprob holds the mixture probabilities for the auxmix
  NumericVector mixprob(10 * T * mpr);
  //mixprob.attr("dim") = Dimension(10, T, mpr, 1); no 4-dim possible?

  // para holds the parameter draws (mu, phi, sigma)
  NumericVector para(3 * mpr * (1));
  para.attr("dim") = Dimension(3, mpr, 1) ;

  // curynorm will hold log((y - facload %*% f)^2) in STEP 1
  arma::mat armaynorm(m, T);

  // curynorm2 will hold log(f^2) in STEP 1
  arma::mat armafnorm(r, T, arma::fill::zeros);
  arma::mat armaXt(r, T);
  arma::mat armaXt2(r, m);
  arma::colvec armaytilde(T);
  arma::colvec armaytilde2(m);
  arma::mat armaSigma(r, r);
  arma::mat armaR(r, r);
  arma::mat armaRinv(r, r);
  arma::mat armaSigma2(r, r);
  arma::mat armaR2(r, r);
  arma::mat armaR2inv(r, r);
  arma::colvec armamean(r);
  arma::colvec armamean2(r);
  arma::colvec armadraw(r);
  arma::colvec armadraw2(r*T);

  // we always use the centered parameterization as baseline
  // (for compatibility reasons with stochvol)
  const bool centered_baseline = true;

  // RNGScope scope;
  // variables are declared afterwards




  // temporary variables for the updated stochvol code
  arma::mat curpara_arma(curpara.begin(), curpara.nrow(), curpara.ncol(), false);
  arma::mat curmixprob_arma(curmixprob.begin(), 10*T, mpr, false);
  arma::imat curmixind_arma(curmixind.begin(), curmixind.nrow(), curmixind.ncol(), false);

  // "linearized residuals"
  // NOTE: "log", "square" are component-wise functions, '*' denotes matrix multiplication
  double offset = 0.00001;
  if (r > 0) {
    armaynorm = log(square(armay - armafacload * armaf));
  } else {
    armaynorm = log(square(armay) + offset);
  }
  armafnorm = log(square(armaf));

  armahtilde = exp(-armah(arma::span::all, arma::span(0,m-1))/2.);


  // STEP 1:
  // update indicators, latent volatilities, and SV-parameters



  // STEP 1 for "linearized residuals"
  for (int j = 0; j < m; j++) {
    if (sv(j) == true) {
      double curh0j = armah0(j);
      arma::vec curpara_j = curpara_arma.unsafe_col(j);
      arma::vec curh_j = armah.unsafe_col(j);
      arma::vec curmixprob_j = curmixprob_arma.unsafe_col(j);
      arma::ivec curmixind_j = curmixind_arma.unsafe_col(j);
      if (j == 0) {
        Rcpp::Rcout << "curmixprob_j: " << curmixprob_j.rows(0, 9).t() << std::endl;
        Rcpp::Rcout << "curmixind_j: " << curmixind_j.rows(0, 9).t() << std::endl;
        Rcpp::Rcout << "curh_j: " << curh_j.rows(0, 9).t() << std::endl;
        Rcpp::Rcout << "data - h: " << arma::mean(armaynorm.row(j)-curh_j.t()) << std::endl;
      }

      stochvol::update_sv(armaynorm.row(j).t(), curpara_j, curh_j, curh0j, curmixprob_j, curmixind_j,
                          centered_baseline, C0(j), cT, Bsigma(j), a0idi, b0idi, bmu, Bmu, B011inv, B022inv, Gammaprior,
                          truncnormal, MHcontrol, MHsteps, parameterization, false, priorh0(j));
      if (j == 0) {
        Rcpp::Rcout << "curmixprob_j after: " << curmixprob_j.rows(0, 9).t() << std::endl;
        Rcpp::Rcout << "curmixind_j after: " << curmixind_j.rows(0, 9).t() << std::endl;
        Rcpp::Rcout << "curh_j after: " << curh_j.rows(0, 9).t() << std::endl;
      }
      armah0(j) = curh0j;

    } else {
      double tmp = sum(square(armay.row(j) - armafacload.row(j)*armaf));
      tmp = 1/as<double>(rgamma(1, priorhomoskedastic(0) + .5*T, 1/(priorhomoskedastic(1) + .5*tmp)));
      armah.col(j).fill(log(tmp));
    }
  }


  // STEP 1 for factors
  for (int j = m; j < mpr; j++) {
    if (sv(j) == true) {
      double curh0j = armah0(j);
      arma::vec curpara_j = curpara_arma.unsafe_col(j);
      arma::vec curh_j = armah.unsafe_col(j);
      arma::vec curmixprob_j = curmixprob_arma.unsafe_col(j);
      arma::ivec curmixind_j = curmixind_arma.unsafe_col(j);
      stochvol::update_sv(armafnorm.row(j-m).t(), curpara_j, curh_j, curh0j, curmixprob_j, curmixind_j,
                          centered_baseline, C0(j), cT, Bsigma(j), a0fac, b0fac, bmu, Bmu, B011inv, B022inv, Gammaprior,
                          truncnormal, MHcontrol, MHsteps, parameterization, true, priorh0(j));
      armah0(j) = curh0j;
    }
  }

  // intermediate step: calculate transformation of curh
  armahtilde = exp(-armah(arma::span::all, arma::span(0,m-1))/2.);

  // STEP 2:
  // update factor loadings: m independent r-variate regressions
  // with T observations (for unrestricted case)
  if (r > 0) {

    int oldpos = 0;
    for (int j = 0; j < m; j++) {

      // TODO: some things outside


      // transposed design matrix Xt is filled "manually"
      int activecols = 0;
      for (int l = 0; l < r; l++) {
          for (int k = 0; k < T; k++) {
            armaXt(activecols, k) = armaf(l, k) * armahtilde(k, j);
          }
          activecols++;
      }

      armaytilde = armay.row(j).t() % armahtilde.col(j);


      // Now draw from the multivariate normal distribution
      // armaSigma is first used as temporary variable:
      armaSigma.submat(0,0,activecols-1,activecols-1) = armaXt.rows(0,activecols-1) * armaXt.rows(0,activecols-1).t();

      // add precisions to diagonal:
      armaSigma.submat(0,0,activecols-1,activecols-1).diag() += 1/arma::nonzeros(armatau2.row(j));

      // Find Cholesky factor of posterior precision matrix
      try {
        armaR.submat(0, 0, activecols-1, activecols-1) = arma::chol(armaSigma.submat(0,0,activecols-1,activecols-1));
      } catch (...) {
        ::Rf_error("Error: Couldn't Cholesky-decompose posterior loadings precision in row %i", j+1);
      }


      // TODO: Check whether Armadillo automatically exploits the fact that R2 is upper triangular for inversion
      // (Partial) Answer: Seems to be OK for native R but solve(trimatu(R), I) is faster with OpenBLAS
      try {
        // armaRinv.submat(0,0,activecols-1,activecols-1) = arma::inv(arma::trimatu(armaR.submat(0,0,activecols-1,activecols-1)));
        armaRinv.submat(0,0,activecols-1,activecols-1) =
          arma::solve(arma::trimatu(armaR.submat(0,0,activecols-1,activecols-1)),
                      arma::eye<arma::mat>(activecols, activecols));
      } catch (...) {
        ::Rf_error("Error: Couldn't invert Cholesky factor of posterior loadings precision in row %i", j+1);
      }

      // calculate posterior covariance armaSigma:
      armaSigma.submat(0, 0, activecols-1, activecols-1) =
        armaRinv.submat(0, 0, activecols-1, activecols-1) *
        armaRinv.submat(0, 0, activecols-1, activecols-1).t();

      // calculate posterior mean:
      armamean.head(activecols) = armaSigma.submat(0, 0, activecols-1, activecols-1) *
        armaXt.submat(0, 0, activecols-1, T-1) *
        armaytilde;

      // draw from the r-variate normal distribution
      armadraw = rnorm(r);

      try {
        armafacloadtmp(arma::span(oldpos, oldpos + activecols - 1)) = armamean.head(activecols) + armaRinv.submat(0,0,activecols-1,activecols-1) * armadraw.head(activecols);
      } catch(...) {
        ::Rf_error("Error: Couldn't sample row %i of factor loadings", j+1);
      }

      //  Rprintf("\n%i to %i: ", oldpos, oldpos+activecols-1);
      //for (int is = oldpos; is < oldpos+activecols; is++) Rprintf("%f ", armafacloadtmp(is));
      //Rprintf("\n\n");
      oldpos = oldpos + activecols;

    }
    armafacloadt(armafacloadtunrestrictedelements) = armafacloadtmp;
    armafacload = arma::trans(armafacloadt);

    //Rprintf("\n\n");
    //for (int is = 0; is < m; is++) Rprintf("%f %f\n", curfacload(is, 0), curfacload(is, 1));
     // STEP 2+: "Deep" Interweaving
      for (int j = 0; j < r; j++) {

        int userow = j;
        if (interweaving == 4) { // find largest absolute element in column to interweave
          userow = 0;
          for (int k = 1; k < m; k++) if (std::fabs(armafacload(k, j)) > std::fabs(armafacload(userow, j))) userow = k;
        }


        //Rprintf("%i and %i\n", j, userow);

        double phi = curpara(1,m+j);
        double sigma = curpara(2,m+j);
        double mu_old = log(armafacload(userow,j) * armafacload(userow,j));
        hopen = armah.col(m+j) + mu_old;
        double h0open = armah0(m+j) + mu_old;
        double logacceptrate;
        double mu_prop;

        if (priorh0(m+j) < 0.) {  // old prior for h0 (stationary distribution, depends on phi), as in JCGS submission Feb 2016
          double tmph = hopen(0) - phi*h0open;
          for (int k = 1; k < T; k++) tmph += hopen(k) - phi*hopen(k-1);

          double gamma_old = (1 - phi) * mu_old;
          double gamma_prop = as<double>(rnorm(1, tmph/(T+B011inv), sigma/std::sqrt(T+B011inv)));
          mu_prop = gamma_prop/(1-phi);

          logacceptrate = logdnormquot(mu_prop, mu_old, h0open, sigma/std::sqrt(1-phi*phi));
          logacceptrate += logspecialquot(gamma_prop, gamma_old, .5, 1/(2.*armatau2(userow,j)), 1-phi);
          logacceptrate += logdnormquot(gamma_old, gamma_prop, 0., sigma*std::sqrt(1/B011inv));

        } else {  // new prior does not depend on phi
          double tmph = hopen(0);
          for (int k = 1; k < (T-1); k++) tmph += hopen(k);

          double tmp4prop = T*priorh0(m+j)*(1-phi)*(1-phi) + 1;
          double prop_mean = (priorh0(m+j) * (1-phi) * (hopen(T-1) + (1-phi)*tmph - phi*h0open) + h0open) / tmp4prop;
          double prop_sd = (sqrt(priorh0(m+j)) * sigma) / std::sqrt(tmp4prop);

          mu_prop = as<double>(rnorm(1, prop_mean, prop_sd));
          logacceptrate = .5 * ((mu_prop - mu_old) - (std::exp(mu_prop) - std::exp(mu_old)) / armatau2(userow,j));
        }

        // NEW, same for both priors:
        arma::vec relevantload = armafacload.col(j);
        arma::vec relevanttau2 = armatau2.col(j);

        // use all except interwoven element (restricted loadings are assumed to be zero!)
        double mysum = accu(square(nonzeros(relevantload))/nonzeros(relevanttau2)) -
          (relevantload(userow)*relevantload(userow))/relevanttau2(userow);

        logacceptrate += .5 * ((nonzerospercol(j)-1)*(mu_prop - mu_old) -
          mysum / (armafacload(userow,j)*armafacload(userow,j)) * (exp(mu_prop) - exp(mu_old)));


        // Rprintf("ACCEPT? ");

        //ACCEPT/REJECT
        if (log(as<double>(runif(1))) < logacceptrate) {
          //    Rprintf("ACC col %i el %02i - ", j+1, userow+1);
          armah.col(m+j) = hopen - mu_prop;
          armah0(m+j) = h0open - mu_prop;

          double tmp = std::exp(mu_prop/2)/armafacload(userow,j);
          armafacload.col(j) *= tmp;
          armaf.row(j) *= 1/tmp;
          //    } else {
          //     Rprintf("REJ col %i el %02i - ", j+1, userow+1);
        }
    }
    // STEP 3:
    // update the factors (T independent r-variate regressions with m observations)

    Rcpp::Rcout << Rcpp::rnorm(1) << std::endl;
    if (samplefac) {
      armadraw2 = rnorm(r*T);
      for (int j = 0; j < T; j++) {

        // transposed design matrix Xt2 (r x m) is filled "manually"
        for (int k = 0; k < m; k++) {
          for (int l = 0; l < r; l++) {
            armaXt2(l, k) = armafacload(k, l) * armahtilde(j,k);
          }
        }

        armaytilde2 = armay.col(j) % armahtilde.row(j).t();

        // Now draw form the multivariate normal distribution

        // armaSigma2 is first used as temporary variable (to hold the precision):
        armaSigma2 = armaXt2 * armaXt2.t();

        // add precisions to diagonal:
        armaSigma2.diag() += exp(-armah(j, arma::span(m, mpr-1)));

        // find Cholesky factor of posterior precision
        try {
          armaR2 = arma::chol(armaSigma2);
        } catch (...) {
          ::Rf_error("Error: Couldn't Cholesky-decompose posterior factor precision at time %i of %i", j+1, T);
        }

        try {
          //   armaR2inv = arma::inv(R2); # This is a little bit faster for very small matrices but a lot slower for large ones...
          //   armaR2inv = arma::inv(arma::trimatu(armaR2)); # This is OK on Native R but not so nice in OpenBLAS
          armaR2inv = arma::solve(arma::trimatu(armaR2), arma::eye<arma::mat>(r, r));
        } catch (...) {
          ::Rf_error("Error: Couldn't invert Cholesky factor of posterior factor precision at time %i of %i",j+1, T);
        }

        // calculate posterior covariance matrix armaSigma2:
        armaSigma2 = armaR2inv * armaR2inv.t();

        // calculate posterior mean armamean2:
        armamean2 = armaSigma2 * armaXt2 * armaytilde2;

        // draw from the r-variate normal distribution
        try {
          armaf.col(j) = armamean2 + (armaR2inv * armadraw2.subvec(j*r, (j+1)*r - 1));
        } catch(...) {
          ::Rf_error("Error: Couldn't sample factors at time %i of %i", j+1, T);
        }
      }
    }
  }


  Rcpp::Rcout << Rcpp::rnorm(1) << std::endl;
  // SIGN SWITCH:
  if (signswitch) {
    for (int j = 0; j < r; j++) {
      if (as<double>(runif(1)) > .5) {
        armafacload.col(j) *= -1;
        armaf.row(j) *= -1;
      }
    }
  }
  Rcpp::Rcout << Rcpp::rnorm(1) << std::endl;
}
