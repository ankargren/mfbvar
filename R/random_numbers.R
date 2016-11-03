rmatn <- function(M, Q, P) {
  L <- chol(Q)
  C <- t(chol(P))
  X <- M + C %*% matrix(rnorm(length(M)), dim(M)) %*% L
  return(X)
}

rmultn <- function(A1, A2) {
  n <- dim(A2)[1]
  X <- A1 + t(chol(A2)) %*% rnorm(n)
  return(X)
}

