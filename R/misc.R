fill_na <- function(Y) {
  apply(Y, 2, function(x) {
    n_x <- length(x) # save lentgh
    if (any(is.na(x))) {
      x <- x[1:max(which(is.na(x) == FALSE))] # get rid of NAs in the end
      for (i in which(is.na(x))) {
        x1 <- NA
        counter <- 1
        while (is.na(x1) == TRUE) {
          x1 <- x[i + counter]
          counter <- counter + 1
        }
        x[i] <- x1
      }

      trimmed_length <- length(x)
      if (trimmed_length < n_x) {
        x <- c(x, rep(NA, n_x - trimmed_length))
        for (i in trimmed_length:n_x) {
          x[i] <- x[trimmed_length]
        }
      }
    }
    x})
}

dnorminvwish <- function(X, Sigma, M, P, S, v) {
  q <- dim(Sigma)[1]
  p <- dim(P)[1]
  det_Sigma <- det(Sigma)
  inv_Sigma <- solve(Sigma)
  dmultnorm <- (-p*q/2) * log(2 * pi) + (-p/2) * log(det_Sigma) + (-q/2)*log(det(P)) + (-1/2 * sum(diag(inv_Sigma %*% t(X - M) %*% solve(P) %*% (X - M))))
  cc <- (v * q/2)*log(2) + (q*(q-1)/4)*log(pi) + sum(lgamma((v+1-1:4)/2))
  dinvwish <- -cc + (v/2) * log(det(S)) -(v+q+1)/2*log(det_Sigma) -1/2 * sum(diag(inv_Sigma %*% S))
  return(exp(dmultnorm + dinvwish))
}

dmultn <- function(x, m, Sigma) {
  p <- dim(Sigma)[1]
  log_d <- (-1/2)* log(det(2*pi*Sigma)) -1/2 * t(x-m) %*% solve(Sigma) %*% (x-m)
  return(exp(log_d))
}

dnorm_trunc <- function(x, m, V_inv, d, p_trunc, chisq_val) {
  qf <- t(x - m) %*% V_inv %*% (x - m)
  return((1/p_trunc) * (1/sqrt((2*pi)^d/det(V_inv))) * exp(-0.5 * qf) * (qf < chisq_val))
}
