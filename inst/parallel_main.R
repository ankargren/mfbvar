
####################################################################
####################################################################

library(parallel)
library(snow)
library(lubridate)
library(mfbvar)
source(system.file("parallel_helpers.R", package = "mfbvar"))
data(mf_list)
start_period <- 49
mf_list$fcst_date <- mf_list$fcst_date[start_period:192]
mf_list$data <- mf_list$data[start_period:192]
mf_list$data <- lapply(mf_list$data, FUN = function(x) x[, -6])
####################################################################
## SETTINGS


lags <- c(4, 12)
freq <- c("QF")
prior_Pi_AR1 <- rep(0, 5)
lambda1_grid <- c(0.1, 0.2)
lambda2_grid <- 1:2
lambda3 <- 10000
n_fcst <- 8
n_burnin <- 100
n_reps <- 100
n_cores <- 1
prior_psi_mean <- c(7, 0.5/3, 0.5, 0, 0.5)
prior_psi_Omega <- diag(0.1, 5)
intervals <- matrix(c(6.5, 7.5,
                      0.4/3, 0.6/3,
                      0.25, 0.75,
                      -1, 1,
                      0.5, 0.65), ncol = 2, byrow = TRUE)
prior_psi_mean <- interval_to_moments(intervals)$prior_psi_mean
prior_psi_Omega <- interval_to_moments(intervals)$prior_psi_Omega

same_seed <- TRUE
save_files <- "C:/Users/seban876/Box Sync/Research/Mixed-Frequency BVAR/results"
cluster_type <- "PSOCK"
data_MF <- mf_list
data_QF <- list(fcst_date = mf_list$fcst_date, data = lapply(mf_list$data, na.omit))

pars <- list(prior_Pi_AR1 = prior_Pi_AR1, lambda1_grid = lambda1_grid, lambda2_grid = lambda2_grid,
                 n_fcst = n_fcst, n_burnin = n_burnin, n_reps = n_reps, prior_psi_mean = prior_psi_mean,
                 prior_psi_Omega = prior_psi_Omega, prior_nu = NULL, lambda3 = lambda3)


####################################################################

####################################################################
## EXPERIMENT 1: MF-SS-4
####################################################################

####################################################################
## EXPERIMENT 2: MF-SS-12
####################################################################

####################################################################
## EXPERIMENT 3: QF-SS-4
####################################################################

####################################################################
## EXPERIMENT 4: QF-SS-12
####################################################################




ex_count <- 0
for (prefix in freq) {
  for (n_lags in lags) {
    for (prior in c("Minn")) {
      ex_count <- ex_count + 1
      Ex_name <- paste("Ex", prefix, n_lags, prior, sep = "_")
      cat("Starting Experiment", Ex_name, paste0("(#", ex_count, ")"), "at", as.character(Sys.time()), "\n")
      if (prefix == "QF") {
        data_list <- data_QF
        pars$Lambda <- build_Lambda(c(rep("identity", 5)), n_lags)
      } else {
        data_list <- data_MF
        pars$Lambda <- build_Lambda(c(rep("identity", 4), "average"), n_lags)
      }

      pars$n_lags <- n_lags

      seed <- 10847 + ex_count

      assign(Ex_name,
             value = parallel_wrapper(data_list, prior, prefix, pars, n_cores, cluster_type = "PSOCK", seed, same_seed))
      #save(list = (Ex_name), file = paste0(save_files, "/", Ex_name, ".RData"))
    }
  }
}

####################################################################
## EXPERIMENT 1b: MF-SS-4 w/ Minnesota
####################################################################

####################################################################
## EXPERIMENT 2b: MF-SS-12 w/ Minnesota
####################################################################

####################################################################
## EXPERIMENT 3b: QF-SS-4 w/ Minnesota
####################################################################

####################################################################
## EXPERIMENT 4b: QF-SS-12 w/ Minnesota
####################################################################

lambda3 <- 1

for (prefix in freq) {
  for (n_lags in lags) {
    ex_count <- ex_count + 1
    cat("Starting Experiment", paste0(ex_count, "b"), "at", as.character(Sys.time()), "\n")
    if (prefix == "QF") {
      data_list <- data_QF
      Lambda <- build_Lambda(c(rep("identity", 6)), n_lags)
    } else {
      data_list <- data_MF
      Lambda <- build_Lambda(c(rep("identity", 4), "average", "identity"), n_lags)
    }

    seed <- 10847 + ex_count

    Ex_num <- paste0("Ex", ex_count, "b")
    assign(Ex_num,
           value = parallel_wrapper_schorf(data_list, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, lambda3,
                                    n_lags, n_fcst, n_burnin, n_reps, n_cores, cluster_type, seed, same_seed))
    save(list = (Ex_num), file = paste0(save_files, "/", Ex_num, ".RData"))
  }
}


