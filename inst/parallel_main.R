if(!("mfbvar" %in% installed.packages()[, "Package"])) {
  if (!("devtools" %in% installed.packages()[, "Package"])) {
    install.packages("devtools")
  }
  devtools::install_github("ankargren/mfbvar")
}

library(mfbvar)
data(mf_list)
data_list <- list()
data_list$data <- lapply(mf_list$data, function(x) {
  x <- as.data.frame(x)
  rownames(x) <- x$date
  x <- x[, -1]
})
data_list$data <- data_list$data[1:6]

n_lags <- 4
Lambda <- build_Lambda(c(rep("identity", 4), "average", "identity"), n_lags)
prior_Pi_AR1 <- rep(0, 6)
lambda1_grid <- c(0.1, 0.2)
lambda2_grid <- 1:2
prior_psi_mean <- c(7, 0.5/3, 1, 0, 0.5, 4)
prior_psi_Omega <- diag(0.1, 6)
n_fcst <- 8
n_burnin <- 100
n_reps <- 100
n_cores <- 6
seed <- 10847

test <- parallel_wrapper(data_list, Lambda, prior_Pi_AR1, lambda1_grid, lambda2_grid, prior_psi_mean, prior_psi_Omega, n_lags, n_fcst, n_burnin, n_reps, n_cores, seed)
