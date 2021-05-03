#' MCMC sampler
#'
#' \code{mcmc_sampler} is a generic function for deciding which specific MCMC
#' algorithm to dispatch to. It is called internally.
#'
#' @param x argument to dispatch on (of class \code{prior_obj})
#' @param ... additional named arguments passed on to the methods
#' @noRd
mcmc_sampler <- function(x, ...) {
  UseMethod("mcmc_sampler")
}

mfbvar_sampler <- function(x, minn = FALSE, ssng = FALSE, ss = FALSE, dl = FALSE,
                       csv = FALSE, fsv = FALSE, iw = FALSE, diffuse = FALSE,
                       ...) {

  envir <- environment()

  params_info <- mfbvar:::get_params_info(minn = minn, ssng = ssng, ss = ss, dl = dl,
                                 csv = csv, fsv = fsv, iw = iw,
                                 diffuse = diffuse)
  mfbvar:::list_to_variables(params_info, envir, "required_params",
                             "prior_params", "retrieved_params", "params")

  ##############################################################################
  ## INPUT CHECK
  mfbvar:::check_required_params(x, required_params)

  ##############################################################################
  ## GET INFO FROM PRIOR OBJECT

  # Assign variables from prior object
  mfbvar:::list_to_variables(x, envir, prior_params)
  cat(prior_params[!(prior_params %in% ls())], "\n")

  # Retrieve some additional variables
  init_vars <- mfbvar:::variable_initialization(Y = Y, freq = freq, freqs = freqs,
                                      n_lags = n_lags, Lambda_ = Lambda_,
                                      n_thin = n_thin)
  mfbvar:::list_to_variables(init_vars, envir, retrieved_params)

  ##############################################################################
  ## PRIOR-SPECIFIC INITIALIZATIONS

  if (iw || csv) {
    priors <- mfbvar:::create_prior_Pi(lambda1 = lambda1,
                                       lambda2 = NULL,
                                       lambda3 = lambda3,
                                       lambda4 = ifelse(minn, lambda4, NULL),
                                       prior_Pi_AR1 = prior_Pi_AR1,
                                       Y = Y,
                                       n_lags = n_lags,
                                       intercept = minn,
                                       prior_nu = n_vars + 2,
                                       independent = FALSE)
    mfbvar:::list_to_variables(priors, envir, "prior_Pi_mean", "prior_Pi_Omega",
                               "prior_S", "inv_prior_Pi_Omega", "Omega_Pi",
                               "prior_nu")
    prior_params <- sort(c(prior_params,
                      "prior_Pi_mean", "prior_Pi_Omega", "prior_S", "prior_nu"))
  }
  if (fsv || diffuse) {
    priors <- mfbvar:::create_prior_Pi(lambda1 = lambda1,
                                               lambda2 = lambda2,
                                               lambda3 = lambda3,
                                               lambda4 = lambda4,
                                               prior_Pi_AR1 = prior_Pi_AR1,
                                               Y = Y,
                                               n_lags = n_lags,
                                               intercept = minn,
                                               block_exo = block_exo,
                                               independent = TRUE)
    mfbvar:::list_to_variables(priors, envir, "prior_Pi_mean", "prior_Pi_Omega")
    prior_params <- sort(c(prior_params,
                           "prior_Pi_mean", "prior_Pi_Omega"))
  }

  # Initalize csv priors
  if (csv) {
    init_csv <- mfbvar:::csv_initialization(prior_phi = prior_phi,
                                            prior_sigma2 = prior_sigma2)
    mfbvar:::list_to_variables(init_csv, envir, "phi_invvar", "phi_meaninvvar",
                               "prior_sigma2", "prior_df", "n_sv")
  }

  # Initialize fsv priors
  if (fsv) {
    init_fsv <- mfbvar:::fsv_initialization(priorsigmaidi = priorsigmaidi,
                                            priorsigmafac = priorsigmafac, priormu = priormu,
                                            priorfacload = priorfacload,
                                            restrict = restrict, priorphiidi = priorphiidi,
                                            priorphifac = priorphifac, n_vars = n_vars,
                                            n_fac = n_fac)
    mfbvar:::list_to_variables(init_fsv, envir, "priorsigmaidi", "priorsigmafac",
                               "bmu", "Bmu", "Bsigma", "B011inv", "B022inv", "armatau2",
                               "armarestr", "a0idi", "b0idi", "a0fac", "b0fac", "priorh0",
                               "n_sv")

  }

  # Initialize ss priors
  if (ss || ssng) {
    init_ss <- mfbvar:::ss_initialization(d = d, d_fcst = d_fcst,
                                          n_T = n_T, n_lags = n_lags, n_fcst = n_fcst)
    mfbvar:::list_to_variables(init_ss, envir, "d_fcst_lags", "D_mat", "dt",
                               "d1", "n_determ")
  }

  # Initialize ssng priors
  if (ssng) {
    init_ssng <- mfbvar:::ssng_initialization(prior_ng, s)
    mfbvar:::list_to_variables(init_ssng, envir, "c0", "c1", "s")
  } else if (ss) {
    phi_mu <- array(0, dim = c(1, 1, 1))
    lambda_mu <- array(0, dim = c(1, 1, 1))
    omega <- array(diag(prior_psi_Omega), dim = c(n_vars * n_determ, 1, n_reps/n_thin))
    c0 <- 0
    c1 <- 0
    s <- 0
    fixate_phi_mu <- FALSE
    fixate_lambda_mu <- FALSE
    fixate_omega <- FALSE
  }

  # Initialize dl priors
  if (dl) {
    init_dl <- dl_initialization(a = a, gig = gig, n_cores = n_cores)
    mfbvar:::list_to_variables(init_dl, envir, "a", "slice", "gig")
  } else if (diffuse || fsv) {
    aux <- array(0, dim = c(1, 1, 1))
    global <- array(0, dim = c(1, 1, 1))
    local <- array(0, dim = c(1, 1, 1))
    a <- -1
    slice <- array(0, dim = c(1, 1, 1))
    gig <- TRUE
    fixate_local <- FALSE
    fixate_aux <- FALSE
    fixate_global <- FALSE
  }

  ##############################################################################
  ## STARTING VALUES FOR PARAMETERS
  add_args <- tryCatch(list(...), error = function(cond) list())
  init <- add_args$init
  fixate <- add_args$fixate
  init_params <- mfbvar:::parameter_initialization(Y = Y, n_vars = n_vars,
                                                   n_lags = n_lags, n_T_ = n_T_,
                                                   init = init, n_fac = n_fac,
                                                   n_determ = n_determ,
                                                   n_sv = n_sv, fsv = fsv,
                                                   csv = csv, params)
  init_fixate <- mfbvar:::fixate_initialization(fixate, params)
  mfbvar:::list_to_variables(init_fixate, envir, paste0("fixate_", params))

  ##############################################################################
  ## INITIALIZATION OF STORAGE OBJECTS
  mfbvar:::storage_initialization(init_params = init_params, params = params,
                                  envir = envir, n_vars = n_vars,
                                  n_lags = n_lags, n_reps = n_reps,
                                  n_thin = n_thin, n_T = n_T, n_T_ = n_T_,
                                  n_determ = n_determ, n_fac = n_fac,
                                  n_fcst = n_fcst, n_sv = n_sv)

  roots <- vector("numeric", n_reps/n_thin)
  num_tries <- roots

  ##############################################################################
  ## SAMPLERS

  # ss(ng) fsv
  if ((ss || ssng) && fsv) {
    mfbvar:::mcmc_ssng_fsv(Y[-(1:n_lags),],Pi,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,
                           mu,phi,sigma,f,facload,latent,
                           Lambda_,prior_Pi_Omega,prior_Pi_AR1,D_mat,dt,d1,
                           d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,bmu,Bmu,
                           a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                           armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,
                           n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,ssng,
                           fixate_Z, fixate_Pi, fixate_psi, fixate_phi_mu,
                           fixate_lambda_mu, fixate_omega, fixate_mu,
                           fixate_phi, fixate_sigma, fixate_f, fixate_facload,
                           fixate_latent)
  }

  # minn/dl fsv
  if (minn && fsv) {
    mfbvar:::mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,latent,
                           aux,global,local,slice,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                           a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                           armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,
                           n_vars,n_T_,n_fcst,n_thin,verbose,a,gig,fixate_Z,fixate_Pi,
                           fixate_mu, fixate_phi, fixate_sigma,
                           fixate_f, fixate_facload, fixate_latent,
                           fixate_aux, fixate_global, fixate_local)
  }

  # ss(ng) iw
  if ((ss || ssng) && iw) {
    mfbvar:::mcmc_ssng_iw(Y[-(1:n_lags),],Pi,Sigma,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,Omega_Pi,prior_Pi_mean,
            prior_S,D_mat,dt,d1,d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,n_reps,n_burnin,
            n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,prior_nu,ssng,
            fixate_Pi, fixate_Sigma, fixate_Z, fixate_psi, fixate_phi_mu,
            fixate_lambda_mu, fixate_omega)
  }

  # minn iw
  if (minn && iw) {
    mfbvar:::mcmc_minn_iw(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,
            Omega_Pi,prior_Pi_mean,prior_S,check_roots,Z_1,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,
            n_thin,verbose,2, fixate_Pi, fixate_Sigma, fixate_Z)
  }

  # ss(ng) csv
  if ((ss || ssng) && csv) {
    mfbvar:::mcmc_ssng_csv(Y[-(1:n_lags),],Pi,Sigma,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,phi,sigma,latent,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,Omega_Pi,prior_Pi_mean,
                  prior_S,D_mat,dt,d1,d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,
                  10,phi_invvar,phi_meaninvvar,prior_sigma2,prior_df,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,ssng, fixate_Pi, fixate_Sigma,
                  fixate_Z, fixate_psi, fixate_phi_mu,
                  fixate_lambda_mu, fixate_omega, fixate_latent,
                  fixate_latent0, fixate_phi, fixate_sigma)
  }

  # minn csv
  if (minn && csv) {
    mfbvar:::mcmc_minn_csv(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,phi,sigma,latent,Lambda_,prior_Pi_Omega,inv_prior_Pi_Omega,
                  Omega_Pi,prior_Pi_mean,prior_S,check_roots,Z_1,10,phi_invvar,phi_meaninvvar,prior_sigma2,prior_df,
                  n_reps,n_burnin,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,n_thin,verbose,fixate_Pi,fixate_Sigma,
                  fixate_Z,fixate_latent,fixate_latent0,fixate_phi,fixate_sigma)
  }

  # minn diffuse
  if (minn && diffuse) {
    mfbvar:::mcmc_minn_diffuse(Y[-(1:n_lags),],Pi,Sigma,Z,Z_fcst,aux,global,local,slice,Lambda_,prior_Pi_Omega,
                      c(prior_Pi_mean),check_roots,Z_1,n_reps,n_burnin,n_q,T_b-n_lags,n_lags,n_vars,n_T_,n_fcst,
                      n_thin,verbose,a,gig,fixate_Pi, fixate_Sigma, fixate_Z,
                      fixate_aux, fixate_global, fixate_local)
  }
  # minn diffuse
  if ((ss || ssng) && diffuse) {
    mfbvar:::mcmc_ssng_diffuse(Y[-(1:n_lags),], Pi, Sigma, psi, phi_mu, lambda_mu,
                               omega, Z, Z_fcst, Lambda_, prior_Pi_Omega, c(prior_Pi_mean), D_mat, dt, d1, d_fcst_lags, prior_psi_mean, c0, c1, s, check_roots, Z_1, n_reps, n_burnin, n_q, T_b-n_lags, n_lags, n_vars, n_T_, n_fcst, n_determ, n_thin, verbose, ssng, fixate_Pi, fixate_Sigma, fixate_Z, fixate_psi, fixate_phi_mu, fixate_lambda_mu, fixate_omega)

  }
  if (verbose) {
    cat("\n")
  }

  ##############################################################################
  ## RETURN OBJECTS AND POST-PROCESSING

  return_obj <- mget(c(params, prior_params, retrieved_params))
  return_obj$add_args <- add_args

  if (check_roots == TRUE) {
    return_obj$roots <- roots
    return_obj$num_tries <- num_tries
  }

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)
}
