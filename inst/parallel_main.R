
####################################################################
####################################################################

library(parallel)
library(snow)
library(lubridate)
library(tidyverse)
library(mfbvar)
source(system.file("parallel_helpers.R", package = "mfbvar"))
data(mf_list)
start_period <- 49
end_period <- 49
mf_list$fcst_date <- mf_list$fcst_date[start_period:end_period]
mf_list$data <- mf_list$data[start_period:end_period]
mf_list$data <- lapply(mf_list$data, FUN = function(x) x[, -6])
####################################################################
## SETTINGS


lags <- c(4)
freq <- c("MF")
prior_Pi_AR1 <- rep(0, 5)
lambda1_grid <- seq(0.05, 1, length.out = 8)
lambda2_grid <- seq(0.5, 5, length.out = 8)
lambda3 <- 10000
n_fcst <- 8
n_burnin <- 5000
n_reps <- 10000
n_cores <- 1
intervals <- matrix(c(6.5, 7.5,
                      0.4/3, 0.6/3,
                      0, 1,
                      -0.1, 0.1,
                      0.5, 0.65), ncol = 2, byrow = TRUE)
prior_psi_mean <- interval_to_moments(intervals)$prior_psi_mean
prior_psi_Omega <- interval_to_moments(intervals)$prior_psi_Omega

same_seed <- TRUE
save_files <- "C:/Users/seban876/Box Sync/Research/Mixed-Frequency BVAR/results"
cluster_type <- "PSOCK"
data_MF <- mf_list
qf_data <- lapply(mf_list$data, FUN = function(x) {
  dat <- x %>%
    mutate(date = rownames(x), year = year(date), quarter = quarter(date)) %>%
    group_by(year, quarter) %>%
    summarize(unemp = mean(unemp), infl = mean(infl), ip = mean(ip), eti = mean(eti), gdp = mean(gdp, na.rm = TRUE)) %>%
    mutate(date = ymd(paste(year, quarter*3, "01", sep = "-")) + months(1)-days(1)) %>%
    ungroup() %>%
    select(-year, -quarter) %>%
    as.data.frame()
  rownames(dat) <- dat$date
  dat %>% select(-date)})
data_QF <- list(fcst_date = mf_list$fcst_date, data = lapply(qf_data, na.omit))

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
    for (prior in c("SS", "Minn")) {
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

Sys.time()
