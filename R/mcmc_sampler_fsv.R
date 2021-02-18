mfbvar_minndl_fsv <- function(x, dl, ...){

  envir <- environment()
  if (dl) {
    required_params <- c("Y", "n_lags", "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                      "n_reps", "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1", "a",
                      "gig", "n_cores")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                          "n_T", "n_T_", "n_thin", "Z_1")
    params <- c("Z", "Pi", "mu", "sigma", "phi", "facload", "f", "latent",
                "latent0", "global", "aux", "local", "slice")

  } else {
    required_params <- c("Y", "n_lags", "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "check_roots", "n_lags", "n_fcst",
                      "n_reps", "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                          "n_T", "n_T_", "n_thin", "Z_1")
    params <- c("Z", "Pi", "mu", "sigma",
                "phi", "facload", "f", "latent", "latent0")
  }


  # Check inputs
  mfbvar:::check_required_params(x, required_params)

  # Assign variables from prior object
  mfbvar:::list_to_variables(x, envir, prior_params)

  # Retrieve some additional variables
  init_vars <- mfbvar:::variable_initialization(Y = Y, freq = freq, freqs = freqs,
                                                n_lags = n_lags, Lambda_ = Lambda_,
                                                n_thin = n_thin)
  mfbvar:::list_to_variables(init_vars, envir, retrieved_params)

  # Prior
  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(lambda1, lambda2, lambda3,
                                                   prior_Pi_AR1, Y, n_lags)

  # Initialize fsv priors
  init_fsv <- mfbvar:::fsv_initialization(priorsigmaidi = priorsigmaidi,
                                          priorsigmafac = priorsigmafac,
                                          priormu = priormu,
                                          priorfacload = priorfacload,
                                          restrict = restrict,
                                          priorphiidi = priorphiidi,
                                          priorphifac = priorphifac,
                                          n_vars = n_vars,
                                          n_fac = n_fac)
  mfbvar:::list_to_variables(init_fsv, envir, "priorsigmaidi", "priorsigmafac",
                             "bmu", "Bmu", "Bsigma", "B011inv", "B022inv",
                             "armatau2", "armarestr", "a0idi", "b0idi", "a0fac",
                             "b0fac", "priorh0")

  # Initialize parameters
  add_args <- tryCatch(list(...), error = function(cond) list())
  init <- add_args$init
  init_params <- mfbvar:::parameter_initialization(Y = Y, n_vars = n_vars, n_lags = n_lags, n_T_ = n_T_,
                                                   init = init, n_fac = n_fac, n_determ = n_determ, params)

  # Initialize storage
  mfbvar:::storage_initialization(init_params = init_params, params = params,
                                  envir = envir, n_vars = n_vars,
                                  n_lags = n_lags, n_reps = n_reps,
                                  n_thin = n_thin, n_T = n_T, n_T_ = n_T_,
                                  n_determ = n_determ, n_fac = n_fac,
                                  n_fcst = n_fcst)

  # Initialize dl
  if (dl) {
    init_dl <- dl_initialization(a = a, gig = gig, n_cores = n_cores)
  } else {
    aux <- matrix(0, 1, 1)
    global <- c(0)
    local <- matrix(0, 1, 1)
    a <- -1
    slice <- c(0)
    gig <- TRUE
  }

  mfbvar:::mcmc_minn_fsv(Y[-(1:n_lags),],Pi,Z,Z_fcst,mu,phi,sigma,f,facload,latent,
                        aux,global,local,slice,Lambda_,prior_Pi_Omega,prior_Pi_AR1, Z_1,bmu,Bmu,
                        a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                        armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags-1,n_lags,
                        n_vars,n_T_,n_fcst,n_thin,verbose,a,gig)
  if (verbose) {
    cat("\n")
  }

  return_obj <- mget(c(params, prior_params, retrieved_params))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }
  return(return_obj)

}

mcmc_sampler.mfbvar_minn_fsv <- function(x, ...){

  out <- mfbvar_minndl_fsv(x, FALSE, ...)

  return(out)
}


mcmc_sampler.mfbvar_dl_fsv <- function(x, ...){

  out <- mfbvar_minndl_fsv(x, TRUE, ...)

  return(out)
}

mfbvar_steadystate_fsv <- function(x, ssng, ...) {

  envir <- environment()
  if (ssng) {
    required_params <- c("Y", "d", "prior_psi_mean", "n_lags",
                          "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "prior_psi_mean",
                      "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                      "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict", "prior_ng", "s",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                      "n_T", "n_T_", "n_thin", "n_determ", "Z_1")
    params <- c("Z", "psi", "Pi", "omega", "phi_mu", "lambda_mu", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")
  } else {
    required_params <- c("Y", "d", "prior_psi_mean", "prior_psi_Omega", "n_lags",
                        "n_burnin", "n_reps", "n_fac")
    prior_params <- c("Y", "freq", "verbose", "prior_psi_mean", "prior_psi_Omega",
                      "d", "d_fcst", "check_roots", "n_lags", "n_fcst", "n_reps",
                      "n_burnin", "n_thin", "n_fac", "freqs", "Lambda_",
                      "priormu", "priorphiidi", "priorphifac", "priorsigmaidi",
                      "priorsigmafac", "priorfacload", "restrict",
                      "lambda1", "lambda2", "lambda3", "prior_Pi_AR1")
    retrieved_params <- c("n_vars", "n_q", "T_b", "n_pseudolags",
                      "n_T", "n_T_", "n_thin", "n_determ", "Z_1")
    params <- c("Z", "psi", "Pi", "mu", "sigma",
              "phi", "facload", "f", "latent", "latent0")
  }


  # Check inputs
  mfbvar:::check_required_params(x, required_params)

  # Assign variables
  mfbvar:::list_to_variables(x, envir, prior_params)

  # Retrieve some additional variables
  init_vars <- mfbvar:::variable_initialization(Y = Y, freq = freq, freqs = freqs,
                                       n_lags = n_lags, Lambda_ = Lambda_,
                                       n_thin = n_thin, d = d, d_fcst = d_fcst)
  mfbvar:::list_to_variables(init_vars, envir, retrieved_params)

  # Prior
  prior_Pi_Omega <- mfbvar:::create_prior_Pi_Omega(lambda1, lambda2, lambda3, prior_Pi_AR1, Y, n_lags)[-1, ]

  # Initialize fsv priors
  init_fsv <- mfbvar:::fsv_initialization(priorsigmaidi = priorsigmaidi,
                                 priorsigmafac = priorsigmafac, priormu = priormu,
                                 priorfacload = priorfacload,
                                 restrict = restrict, priorphiidi = priorphiidi,
                                 priorphifac = priorphifac, n_vars = n_vars,
                                 n_fac = n_fac)
  mfbvar:::list_to_variables(init_fsv, envir, "priorsigmaidi", "priorsigmafac",
                    "bmu", "Bmu", "Bsigma", "B011inv", "B022inv", "armatau2",
                    "armarestr", "a0idi", "b0idi", "a0fac", "b0fac", "priorh0")

  # Initialize ss objects
  init_ss <- mfbvar:::ss_initialization(d = d, d_fcst = d_fcst,
                               n_T = n_T, n_lags = n_lags, n_fcst = n_fcst)
  mfbvar:::list_to_variables(init_ss, envir, "d_fcst_lags", "D_mat", "dt", "d1")

  # Initialize parameters
  add_args <- tryCatch(list(...), error = function(cond) list())
  init <- add_args$init
  init_params <- mfbvar:::parameter_initialization(Y = Y, n_vars = n_vars, n_lags = n_lags, n_T_ = n_T_,
                  init = init, n_fac = n_fac, n_determ = n_determ, params)

  # Initialize storage
  mfbvar:::storage_initialization(init_params = init_params, params = params, envir = envir,
                         n_vars = n_vars, n_lags = n_lags, n_reps = n_reps,
                         n_thin = n_thin, n_T = n_T, n_T_ = n_T_,
                         n_determ = n_determ, n_fac = n_fac, n_fcst = n_fcst)

  # Initialize ssng priors
  if (ssng) {
    init_ssng <- mfbvar:::ssng_initialization(prior_ng, s)
    mfbvar:::list_to_variables(init_ssng, envir, "c0", "c1", "s")
  } else {
    phi_mu <- matrix(0, 1, 1)
    lambda_mu <- matrix(0, 1, 1)
    omega <- matrix(diag(prior_psi_Omega), nrow = 1)
    c0 <- 0
    c1 <- 0
    s <- 0
  }

  roots <- vector("numeric", n_reps/n_thin)
  num_tries <- roots

  mfbvar:::mcmc_ssng_fsv(Y[-(1:n_lags),],Pi,psi,phi_mu,lambda_mu,omega,Z,Z_fcst,
                       mu,phi,sigma,f,facload,latent,
                       Lambda_,prior_Pi_Omega,prior_Pi_AR1,D_mat,dt,d1,
                       d_fcst_lags,prior_psi_mean,c0,c1,s,check_roots,Z_1,bmu,Bmu,
                       a0idi,b0idi,a0fac,b0fac,Bsigma,B011inv,B022inv,priorh0,
                       armarestr,armatau2,n_fac,n_reps,n_burnin,n_q,T_b-n_lags-1,n_lags,
                       n_vars,n_T_,n_fcst,n_determ,n_thin,verbose,ssng)
  if (verbose) {
    cat("\n")
  }

  return_obj <- mget(c(params, prior_params, retrieved_params))

  if (n_fcst > 0) {
    return_obj$Z_fcst <- Z_fcst
  }

  return(return_obj)
}
mcmc_sampler.mfbvar_ss_fsv <- function(x, ...){

  out <- mfbvar_steadystate_fsv(x, FALSE, ...)

  return(out)
}

mcmc_sampler.mfbvar_ssng_fsv <- function(x, ...){

  out <- mfbvar_steadystate_fsv(x, TRUE, ...)

  return(out)
}
