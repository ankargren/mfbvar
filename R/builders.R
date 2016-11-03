# Contains:
# D_f
# Z_tilde_f
# U_f

build_DD <- function(d, n_lags) {
  # Inputs:

  n_determ <- ncol(d)
  n_T <- nrow(d)
  DD <- d[-(1:(n_lags)), ]
  for (i in 1:n_lags) {
    DD <- cbind(DD, -d[-c(0:(n_lags-i), ((n_T-i+1):n_T)), ])
  }
  return(DD)
  # The output is a (t_1-t_0+1)*pk matrix
}


build_demeaned_z <- function(z, psi, d) {
  # Inputs:
  # t_0 (scalar), start period
  # t_1 (scalar), end period
  # Z is a matrix, T * p
  # psi is mp * 1, vectorized Psi

  n_vars <- ncol(z)
  n_T <- nrow(z)
  demeaned_z <- z - kronecker(diag(n_T), matrix(psi, nrow = 1)) %*% kronecker(matrix(c(t(d)), ncol = 1), diag(n_vars))
  return(demeaned_z)
  # The output is a (t_1-t_0+1)*pk matrix
}

build_Z <- function(z, n_lags) {
  # Inputs:

  n_vars <- ncol(z)
  n_T <- nrow(z)
  Z <- z[-(1:(n_lags-1)), ]
  for (i in 2:n_lags) {
    Z <- cbind(Z, z[-c(0:(n_lags-i), ((n_T-i+2):n_T)), ])
  }
  return(Z)
  # The output is a (t_1-t_0+1)*pk matrix
}


build_U <- function(Pi, n_determ) {
  # Pi is (Pi_1, ..., Pi_k)'
  n_vars <- dim(Pi)[1]
  n_lags <- dim(Pi)[2]/n_vars
  U <- diag(n_vars*n_determ)
  for(i in 1:n_lags) {
    U <- rbind(U, kronecker(diag(n_determ), Pi[, (1 + (i-1)*n_vars):(i*n_vars)]))
  }
  return(U)
}

build_Y_tilde <- function(Pi, z) {
  n_vars <- nrow(Pi)
  n_lags <- ncol(Pi)/n_vars
  Z <- build_Z(z, n_lags) # This gives Z_{0:T}, we need Z_{0:T-1}
  Z <- Z[-nrow(Z), ]

  Y_tilde <- z[-(1:(n_lags)), ] - Z %*% t(Pi)
  return(Y_tilde)
}

build_M_Lambda <- function(Y, Lambda, n_lags) {
  M_Lambda <- list()
  n_T <- nrow(Y)
  for (i in 1:n_T) {
    M_Lambda[[i]] <- matrix(diag(n_vars), ncol = n_vars) %*% Lambda
    if (any(is.na(Y[i, ]))) {
      M_Lambda[[i]][is.na(Y[i,]), ] <- NA
    }
  }
  return(M_Lambda)
}



