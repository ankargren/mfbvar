library(Rcpp)
library(RcppArmadillo)
load("y.RData")
y <- cbind(y[, 2], y[, 1])
y <- cbind(y, y + matrix(rnorm(nrow(y)*2), ncol =2))
y <- cbind(y[, c(1, 3)], y[, c(2, 4)])

# Load files
files <- c("builders", "prior_pi_sigma", "posteriors", "ols", "misc", "smooth_samp", "random_numbers")
sapply(files, function(x) source(paste0(getwd(), "/R/", x, ".R")))
sourceCpp(paste0(getwd(), "/src/rnd_numbers.cpp"))
sourceCpp(paste0(getwd(), "/src/builders.cpp"))

# Set priors
priors <- prior_pi_sigma(lambda1 = 0.2, lambda2 = 1, prior_mean = c(0.5, 0.5, 0.5, 0.5), Y = y,
                               n_lags = 4, nu = 4)
prior_pi <- priors$prior_pi
prior_pi_omega <- priors$prior_pi_omega
prior_nu <- 4
prior_s <- priors$prior_s
prior_psi <- c(0.2, 0.4, 0, 0, 0.2, 0.1, 0.05, 0.08)
prior_psi_omega <- diag(0.5, 8)

Y <- y[-c(1:12), ]
d <- matrix(cbind(1, 1:200), dim(Y))

Lambda <- rbind(c(1, rep(0, 7)), c(0, rep(c(1/3, 0), 3), rep(0, 1)))
Lambda <- rbind(Lambda[c(1, 1),], Lambda[c(2, 2),])
n_lags <- 4
n_vars <- 4
M_Lambda <- build_M_Lambda(Y[-(1:n_lags), ], Lambda, n_lags)
n_burnin <- 100
n_reps <- 100



gibbs_sampler <- function(prior_pi, prior_pi_omega, prior_nu, prior_s, prior_psi, prior_psi_omega,
                          Y, d, M_Lambda, n_reps = 10000, n_burnin = 10000) {
  # prior_pi: (p * kp) matrix of prior mean for Pi
  # prior_pi_omega: (kp^2 * kp^2) matrix of prior covariance matrix for vec(pi')
  # prior_nu: scalar with the prior for nu
  # prior_s: (p * p) matrix with the prior for s
  # prior_psi:
  # prior_psi_omega:

  # Y: (T * p) matrix of main data (with NA where observations are missing)
  # d: (T * m) matrix of deterministic data
  # Lambda:


  # n_burnin: scalar with the number of burn-in replications
  # n_lags: scalar with the number of lags
  # n_reps: scalar with the number of replications

  n_vars <- dim(Y)[2]
  n_lags <- prod(dim(as.matrix(prior_pi)))/n_vars^2
  n_determ <- dim(d)[2]
  n_tot_reps <- n_reps + n_burnin
  n_T_full <- dim(Y)[1]
  n_T <- n_T_full# - n_lags

  #################
  # Preallocation
  # Pi and Sigma store their i-th draws in the third dimension, psi
  # is vectorized so it has its i-th draw stored in the i-th row
  # Pi:    p * pk * n_tot_reps, each [,,i] stores Pi'
  # Sigma: p * p  * n_tot_reps
  # psi:   n_tot_reps * p
  # pre_Z: (T_ + k + 1) * p * n_tot_reps
  # Z:     T * p * n_tot_reps
  Pi    <- array(NA, dim = c(n_vars, n_vars * n_lags, n_tot_reps))
  Sigma <- array(NA, dim = c(n_vars, n_vars, n_tot_reps))
  psi   <- array(NA, dim = c(n_tot_reps, n_vars * n_determ))
  Z     <- array(NA, dim = c(n_T, n_vars, n_tot_reps))

  #################
  # Initialization
  # The missing values in Z are filled with the next observed value
  # Pi, Sigma and psi are then computed using maximum likelihood

  Z[,, 1] <- fill_na(Y)

  ols_results <- ols_initialization(Z[,, 1], d = d, n_lags = n_lags)

  Pi[,, 1]    <- ols_results$Gam[, 1:(n_vars * n_lags)]
  Sigma[,, 1] <- ols_results$S
  psi[1, ]    <- c(solve(diag(n_vars) - Pi[,, 1] %*%
                        kronecker(matrix(1, n_lags, 1), diag(n_vars))) %*%
                  ols_results$Gam[, (n_vars * n_lags + 1):(n_vars * n_lags + n_determ)])
  D <- build_DD(d, n_lags)


  for (r in 2:(n_tot_reps)) {
    #demeaned_Z <- build_demeaned_z(Z[,,r-1], psi[r-1, ], d) %>%
      #build_Z(n_lags)

    demeaned_Z <- build_demeaned_z_cpp(Z[,,r-1], matrix(psi[r-1, ], nrow = 1), matrix(d, ncol = n_determ), n_T, n_vars) %>%
      build_Z(n_lags)

    # Dynamic coefficients
    post_pi_omega <- demeaned_Z[-nrow(demeaned_Z), ] %>% posterior_pi_omega(prior_pi_omega, .)
    post_pi   <- posterior_pi(post_pi_omega, prior_pi_omega, prior_pi, demeaned_Z)
    pi_sample <- ols_pi(demeaned_Z[-nrow(demeaned_Z), ], demeaned_Z[-1, 1:n_vars])

    # Covariance
    s_sample  <- ols_s(demeaned_Z[-nrow(demeaned_Z), ], demeaned_Z[-1, 1:n_vars], pi_sample)
    post_s <- posterior_s(prior_s, s_sample, prior_pi, pi_sample, prior_pi_omega, demeaned_Z[-nrow(demeaned_Z), ])
    nu <- n_T + prior_nu # Is this the right T? Or should it be T - lags?
    Sigma[,,r] <- rinvwish(nu, post_s)

    # Pi
    Pi[,,r] <- rmatnorm(post_pi, Sigma[,,r], post_pi_omega)

    # psi
    U <- build_U_cpp(Pi[,,r], n_determ, n_vars, n_lags)
    post_psi_omega <- posterior_psi_omega(U, D, Sigma[,, r], prior_psi_omega)
    Y_tilde <- build_Y_tilde(Pi[,, r], Z[,, r-1])
    post_psi <- posterior_psi(post_psi_omega, U, Sigma[,, r], Y_tilde, D, prior_psi_omega, prior_psi)
    psi[r, ] <- t(rmultn(post_psi, post_psi_omega))

    # Z

    Pi_comp    <- rbind(Pi[,, r], cbind(diag(n_vars*(n_lags-1)), matrix(0, ncol = n_vars, nrow = n_vars*(n_lags-1))))
    psi_comp   <- rbind(matrix(psi[r, ], ncol = n_determ, nrow = n_vars),
                       matrix(0, nrow = n_vars*(n_lags - 1), ncol = n_determ))
    Q_comp     <- matrix(0, ncol = n_vars*n_lags, nrow = n_vars*n_lags)
    Q_comp[1:n_vars, 1:n_vars] <- t(chol(Sigma[,,r]))


    mZ  = Y[-(1:n_lags), ]
    mX  = d[-(1:n_lags), ]
    lH  = M_Lambda
    lH0 = NULL
    mF  = Pi_comp
    mB  = psi_comp
    mQ  = Q_comp
    iT  = n_T - n_lags
    ip  = n_vars
    iq  = n_lags * n_vars
    is  = n_determ
    h0  = matrix(c(t(Z[1:n_lags,,1])), ncol = 1)
    P0  = NULL
    X0  = matrix(d[n_lags, ], ncol = 1)

    smoothed_Z <- smooth_samp(mZ,mX,lH,lH0=NULL,mF,mB,mQ,iT,ip,iq,is,h0,P0=NULL,X0)

    # Forecasting
  }

}





