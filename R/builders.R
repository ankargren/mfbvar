#' Build the \eqn{D} matrix
#'
#' \code{build_DD} builds the \eqn{D} matrix.
#'
#' @templateVar d TRUE
#' @templateVar n_lags TRUE
#' @template man_template
#'
#' @return
#' \item{DD}{A matrix of size \code{n_T * ((n_lags + 1)*n_determ)} where row \code{t} is \eqn{(d_t', -d_{t-1}', \dots, -d_{t-k}')}.}

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

#' Build the companion matrix for the dynamic parameters
#'
#' Builds the parameter matrix of dynamic coefficients for the companion form representation.
#'
#' @templateVar Pi TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUE
#' @template man_template
#'
#' @return
#' \item{Pi_comp}{The companion form matrix of size \code{(n_vars*n_lags) * (n_vars*n_lags)}}.

build_companion <- function(Pi, n_vars, n_lags) {
  rbind(Pi, cbind(diag(n_vars*(n_lags-1)), matrix(0, ncol = n_vars, nrow = n_vars*(n_lags-1))))
}

#' Build the \eqn{Z} matrix
#'
#' Builds the \eqn{Z} matrix, which consists of lags of \eqn{z}.
#'
#' @templateVar z TRUE
#' @templateVar n_lags TRUE
#' @template man_template
#'
#' @return
#' \item{Z}{A matrix of size \code{n_T * (n_vars*n_lags)}.}

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

#' Build the \code{U} matrix
#'
#' Builds the parameter matrix of dynamic coefficients for the companion form representation.
#'
#' @templateVar Pi TRUE
#' @templateVar n_determ TRUE
#' @template man_template
#'
#' @return
#' \item{U}{The \code{U} matrix, of size \code{((n_lags+1)n_vars*n_determ) * n_vars*n_determ}.}

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

#' Build the lag-corrected data matrix
#'
#' Builds the \eqn{\tilde{Y}=\Pi(L)Y} matrix.
#'
#' @templateVar Pi TRUE
#' @templateVar z TRUE
#' @template man_template
#'
#' @return
#' \item{Y_tilde}{A matrix of size \code{n_T * n_vars}.}
#' @details Note that \code{z} does not contain missing values; at this point, the missing values have been replaced by values drawn using the simulation smoother.

build_Y_tilde <- function(Pi, z) {
  n_vars <- nrow(Pi)
  n_lags <- ncol(Pi)/n_vars
  Z <- build_Z(z, n_lags) # This gives Z_{0:T}, we need Z_{0:T-1}
  Z <- Z[-nrow(Z), ]

  Y_tilde <- z[-(1:(n_lags)), ] - Z %*% t(Pi)
  return(Y_tilde)
}

#' Build the \eqn{M_t\Lambda} matrices
#'
#' Builds the selection matrices \eqn{M_t\Lambda}.
#' @templateVar Y TRUE
#' @templateVar Lambda TRUE
#' @templateVar n_vars TRUE
#' @templateVar n_lags TRUE
#' @templateVar n_T TRUE
#' @template man_template
#'
#' @return
#' \item{M_Lambda}{A list of length \code{n_T}.}
#' @details The element \code{M_Lambda[[t]]} corresponds to \eqn{M_t\Lambda}. Currently, if element \code{i} of \code{Y[t, ]} is \code{NA}, then row \code{i} of \code{M_Lambda[[t]]} is all \code{NA}.

build_M_Lambda <- function(Y, Lambda, n_vars, n_lags, n_T) {
  M_Lambda <- list()
  for (i in 1:n_T) {
    M_Lambda[[i]] <- matrix(diag(n_vars), ncol = n_vars) %*% Lambda
    if (any(is.na(Y[i, ]))) {
      M_Lambda[[i]][is.na(Y[i,]), ] <- NA
    }
  }
  return(M_Lambda)
}

#' Build the \eqn{\Lambda} matrix
#'
#' Builds the aggregation matrix \eqn{\Lambda}.
#' @templateVar aggregation TRUE
#' @templateVar n_lags TRUE
#' @template man_template
#' @return
#' \item{Lambda}{An \code{n_vars * (n_vars*n_pseudolags)} matrix, where \code{n_pseudolags} is \code{max(5, n_lags)} if any variable uses the triangular aggregation scheme, \code{max(3, n_lags)} if any uses the quarterly average.}
#' @details The choice \code{aggregation = "identity"} means that what is observed is assumed to be exactly the underlying variable, whereas \code{aggregation = "average"} uses the quarterly average of the monthly underlying observations. Lastly, \code{aggregation = "triangular"} uses the triangular specification used by Mariano and Murasawa (2010).
build_Lambda <- function(aggregation, n_lags) {
  n_vars <- length(aggregation)
  if (any(aggregation %in% "triangular") && n_lags < 5) {
    Lambda <- matrix(0, n_vars, n_vars * 5)
  } else if (n_lags > 2) {
    Lambda <- matrix(0, n_vars, n_vars * n_lags)
  } else {
    stop("Too few lags!")
  }

  n_pseudolags <- dim(Lambda)[2]/n_vars
  for (i in 1:n_vars) {
    if (aggregation[i] == "identity") {
      fill_vec <- c(1, rep(0, n_pseudolags - 1))
    }
    if (aggregation[i] == "average") {
      fill_vec <- c(rep(1/3, 3), rep(0, n_pseudolags - 3))
    }
    if (aggregation[i] == "triangular") {
      fill_vec <- c(1/3, 2/3, 1, 2/3, 1/3, rep(0, n_pseudolags - 5))
    }

    Lambda[i, seq(i, n_pseudolags * n_vars, by = n_vars)] <- fill_vec
  }
  return(Lambda)
}


