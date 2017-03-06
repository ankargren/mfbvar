library(expm)
library(MFBVAR)
Phi <- cbind(matrix(c(0.8, 0.1, 0.2, 0.5), 2, 2), matrix(0, 2, 6))
Sigma <- diag(c(1, 2))
n_vars <- 2
n_lags <- 4
n_T <- 100

linear_transformation <- function() {
  stacked <- c(t(Y))
  identity <- diag(n_vars*n_T)
  E0 <- identity[ is.na(stacked), ]
  E1 <- identity[!is.na(stacked), ]
  
  A <- kronecker(diag(n_T/n_lags), A_tilde)
  B <- kronecker(diag(n_T/n_lags), kronecker(diag(n_lags)[1,,drop=FALSE], diag(n_vars*n_lags)))
  
  Phi_comp <- build_companion(Phi, n_vars = n_vars, n_lags = n_lags)
  Sigma_comp <- matrix(0, n_vars*n_lags, n_vars*n_lags)
  Sigma_comp[1:n_vars, 1:n_vars] <- Sigma
  
  U <- lapply(1:n_T, function(x) (Phi_comp %^% (x-1)) %*% Sigma_comp)
  
  k <- min(unlist(lapply(U, dim)))
  n <- length(U)
  strip <- array(NA, dim=c(k,k,2*n-1))
  for (i in 1:n) strip[,,i] <- U[[n+1-i]]
  if (n > 1) for (i in 2:n) strip[,,n+i-1] <- t(U[[i]])
  X <- array(NA, dim=c(k,k,n,n))
  for (i in 1:n) X[,,,i] <- strip[,,(n+1-i):(2*n-i)]
  Omega <- matrix(aperm(X, c(1,3,2,4)), n*k)
  
  
}