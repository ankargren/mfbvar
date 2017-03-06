library(expm)
library(MFBVAR)

## Generate some data. Two variables, 200 time points.
TT <- 200
set.seed(100)

Y <- matrix(0, 2*TT, 2)
Phi <- matrix(c(0.5, 0.1, 0.2, 0.5), 2, 2)
for (i in 2:(2*TT)) {
  Y[i, ] <- Phi %*% Y[i-1,] + rnorm(2)
}
Y[, 2] <- zoo::rollapply(Y[, 2], 3, mean, fill = NA, align = "right")
Y <- Y[-(1:TT),]
Y[setdiff(1:TT, seq(1, TT, 3)), 2] <- NA
Y <- Y[-c(1, 200), c(2, 1)]


Phi <- cbind(matrix(c(0.8, 0.1, 0.2, 0.5), 2, 2), matrix(0, 2, 6))
Sigma <- diag(c(1, 2))
n_vars <- 2
n_lags <- 4
n_T <- 198
A_tilde <- diag(6)
A_tilde[5,] <- rep(c(1/3, 0), 3)

# This assumes that Y starts at beginning of quarter!
linear_transformation <- function(Y, Phi, Sigma, A_tilde, n_vars, n_T, n_lags) {
  stacked <- c(t(Y))
  identity <- diag(n_vars*n_T)
  E0 <- identity[ is.na(stacked), ]
  E1 <- identity[!is.na(stacked), ]
  
  n_agglags <- (nrow(A_tilde)/n_vars)
  A <- kronecker(diag(ceiling(n_T/n_agglags)), A_tilde)[1:(n_vars*n_T), 1:(n_vars*n_T)]
  # A: Z_1, Z_2, ...
  # A <- A[nrow(A):1, ncol(A):1]
  # B <- kronecker(diag(ceiling((n_T-n_lags+1)/n_lags)), kronecker(diag(n_lags)[1,,drop=FALSE], diag(n_vars*n_lags)))[, 1:((n_T-n_lags+1) * n_lags * n_vars)]#[1:(floor(n_T/n_lags)*n_vars*n_lags + (n_T-floor(n_T/n_lags)*n_lags)*n_vars), ]
  # missing_rows <- n_vars*n_T-dim(B)[1]
  # B <- rbind(B, matrix(0, nrow = missing_rows, ncol = dim(B)[2]))
  # B[(nrow(B)-missing_rows+1):nrow(B), (ncol(B)-missing_rows+1):ncol(B)] <- diag(missing_rows)
  
  B <- kronecker(diag(n_T), kronecker(diag(n_lags)[1,,drop=FALSE], diag(n_vars)))
  # B works with a vector (z_T, z_(t-1), ..., z_p)
  # c(t(build_Z(Y, 4)[8:1,]))
  #build_Z(rbind(matrix(0, 3, 2), Y), 4)
  Phi_comp <- build_companion(Phi, n_vars = n_vars, n_lags = n_lags)
  Sigma_comp <- matrix(0, n_vars*n_lags, n_vars*n_lags)
  Sigma_comp[1:n_vars, 1:n_vars] <- Sigma
  
  U <- lapply(1:(n_T), function(x) (Phi_comp %^% (x-1)) %*% Sigma_comp)
  
  k <- min(unlist(lapply(U, dim)))
  n <- length(U)
  strip <- array(NA, dim=c(k,k,2*n-1))
  for (i in 1:n) strip[,,i] <- U[[n+1-i]]
  if (n > 1) for (i in 2:n) strip[,,n+i-1] <- t(U[[i]])
  X <- array(NA, dim=c(k,k,n,n))
  for (i in 1:n) X[,,,i] <- strip[,,(n+1-i):(2*n-i)]
  Omega <- matrix(aperm(X, c(1,3,2,4)), n*k)
  
  AB <- A %*% B
  E0AB <- E0 %*% AB
  E1AB <- E1 %*% AB
  E0ABOmegatE0AB <- E0AB %*% Omega %*% t(E0AB)
  E0ABOmegatE1AB <- E0AB %*% Omega %*% t(E1AB)
  E1ABOmegatE1AB <- E1AB %*% Omega %*% t(E1AB)
  C <- E0ABOmegatE1AB %*% chol2inv(chol(E1ABOmegatE1AB))
  mu_bar <- C %*% stacked[!is.na(stacked)]
  Omega_bar <- E0ABOmegatE0AB - C %*% t(E0ABOmegatE1AB)
  
  draw <- rmultn(m = mu_bar, Sigma = Omega_bar)
  # draw is x_1, x_2, x_4, x_5, ...
  A_tilde_inv <- solve(A_tilde[n_vars*(1:n_agglags)-1, n_vars*(1:n_agglags)-1])
  A_inv <- kronecker(diag(ceiling(n_T/n_agglags)), A_tilde_inv)[1:n_T, 1:n_T]
  
  Y[is.na(Y[, 1]), 1] <- draw
  Y[, 1] <- A_inv %*% Y[, 1]
  return(Y)
}

temp <- array(Y, dim = c(nrow(Y), ncol(Y), 100))
for (i in 1:100) {
  temp[,,i] <- linear_transformation(Y, Phi, Sigma, A_tilde, n_vars, n_T, n_lags)
}

temp <- temp[,1,]

library(tidyverse)
plot_df <- as_tibble(Y) %>%
  select(Y1 = V1) %>%
  mutate(lower  = apply(temp, 1, quantile, prob = 0.05),
         median = apply(temp, 1, quantile, prob = 0.5),
         upper  = apply(temp, 1, quantile, prob = 0.95),
         time = 1:n()) %>%
  gather(Y1, median, key = "key", value = "value") %>%
  na.omit()
ggplot() +
  geom_ribbon(data = plot_df, aes(x = time, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = filter(plot_df, key %in% c("Y1", "median")), aes(x = time, y = value, color = key))
  

mZ <- Y[-(1:n_lags), ]
Lambda <- matrix(0, 2, 8)
Lambda[1, ] <- c(rep(c(1/3, 0), 3), 0, 0)
Lambda[2, 2] <- 1
mF <- build_companion(Phi, n_vars, n_lags)
mQ <- matrix(0, n_vars*n_lags, n_vars*n_lags)
mQ[1:n_vars, 1:n_vars] <- t(chol(Sigma))
iT <- n_T-n_lags
ip <- n_vars
iq <- n_vars*n_lags
h0 <- Y[1:n_lags, ] 
h0[c(1, 2, 4), 1] <- h0[3, 1]
h0 <- matrix(h0, ncol = 1)
P0 <- mQ * 0


temp2 <- array(Y, dim = c(nrow(Y), ncol(Y), 100))
for (i in 1:100) {
  temp2[,,i] <- rbind(matrix(h0, ncol = 2), simulation_smoother(mZ, Lambda, mF, mQ, iT, ip, iq, h0, P0)[, 1:2])

}

temp2 <- temp2[,1,]

library(tidyverse)
plot_df <- as_tibble(Y) %>%
  select(Y1 = V1) %>%
  mutate(lower  = apply(temp2, 1, quantile, prob = 0.05),
         median = apply(temp2, 1, quantile, prob = 0.5),
         upper  = apply(temp2, 1, quantile, prob = 0.95),
         time = 1:n()) %>%
  gather(Y1, median, key = "key", value = "value") %>%
  na.omit()
ggplot() +
  geom_ribbon(data = plot_df, aes(x = time, ymin = lower, ymax = upper), alpha = 0.3) +
  geom_line(data = filter(plot_df, key %in% c("Y1", "median")), aes(x = time, y = value, color = key))

plot_df <- as_tibble(Y) %>%
  select(Y1 = V1) %>%
  mutate(LF  = apply(temp, 1, quantile, prob = 0.5),
         KF = apply(temp2, 1, quantile, prob = 0.5),
         time = 1:n()) %>%
  gather(Y1, LF, KF, key = "key", value = "value") %>%
  na.omit()
ggplot() +
  geom_line(data = plot_df, aes(x = time, y = value, color = key))

