if(!("mfbvar" %in% installed.packages()[, "Package"])) {
  if (!("devtools" %in% installed.packages()[, "Package"])) {
    install.packages("devtools")
  }
  devtools::install_github("ankargren/mfbvar")
}

list.of.packages <- c("parallel", "lubridate", "snow", "Rmpi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(parallel)
library(snow)
library(Rmpi)
library(lubridate)
library(mfbvar)
source(system.file("parallel_helpers.R", package = "mfbvar"))
data(mf_list)


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
lags <- c(4, 12)
freq <- c("MF", "QF")
prior_Pi_AR1 <- rep(0, 6)
lambda1_grid <- c(0.1, 0.2)
lambda2_grid <- 1:2
n_fcst <- 8
n_burnin <- 100
n_reps <- 100
n_cores <- 6
for (prefix in freq) {
  for (n_lags in lags) {
    ex_count <- ex_count + 1
    cat("Starting Experiment", ex_count, "at", as.character(Sys.time()))
    data_list <- list()
    data_list$data <- mf_list$data[1:6]
    if (prefix == "QF") {
      data_list$data <- lapply(data_list$data, na.omit)
      Lambda <- build_Lambda(c(rep("identity", 6)), n_lags)
      prior_psi_mean <- c(7, 0.5/3, 1, 0, 0.5, 4)
      prior_psi_Omega <- diag(0.1, 6)
    } else {
      Lambda <- build_Lambda(c(rep("identity", 4), "average", "identity"), n_lags)
      prior_psi_mean <- c(7, 0.5/3, 1, 0, 0.5, 4)
      prior_psi_Omega <- diag(0.1, 6)
    }

    seed <- 10847 + ex_count

    Ex_num <- paste0("Ex", ex_count)
    assign(Ex_num,
           value = parallel_wrapper(data_list, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, prior_psi_mean, prior_psi_Omega,
                                    n_lags, n_fcst, n_burnin, n_reps, n_cores, seed))
    save(list = (Ex_num), file = paste0(Ex_num, ".RData"))
  }
}

