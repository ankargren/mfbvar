
list_to_matrix <- function(Y_in) {
  if (all(sapply(Y_in, function(x) inherits(x, "ts"))) || all(sapply(Y_in, function(x) inherits(x, "zoo")))) {
    if (all(sapply(Y_in, function(x) inherits(x, "ts")))) {
      zoofun <- function(x) {
        if (frequency(x) == 4) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date.ts(x) %m+% months(2)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date.ts(x) %m+% months(2)))
          }
        } else if (frequency(x) == 12) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date.ts(x)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date.ts(x)))
          }
        } else {
          stop("Time series objects can only include monthly and/or quarterly time series.")
        }
      }
    } else if (all(sapply(Y_in, function(x) inherits(x, "zooreg")))) {
      zoofun <- function(x) {
        if (frequency(x) == 4) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date(zoo::index(x)) %m+% months(2)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date(zoo::index(x)) %m+% months(2)))
          }
        } else if (frequency(x) == 12) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date(zoo::index(x))))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date(zoo::index(x))))
          }
        } else {
          stop("Time series objects can only include monthly and/or quarterly time series.")
        }
      }
    }
    zoolist <- lapply(Y_in, zoofun)
    reducedlist <- Reduce(zoo::merge.zoo, zoolist)
    Y <- as.matrix(reducedlist)
    rownames(Y) <- as.character(time(reducedlist))
    dim_null <- sapply(zoolist, function(x) is.null(dim(x)))
    if (all(dim_null)) {
      colnames(Y) <- names(zoolist)
    } else if (all(!dim_null)) {
      colnames(Y) <- Reduce(c, lapply(zoolist, colnames))
    } else {
      name_vec <- c()
      for (iter in 1:length(dim_null)) {
        if (dim_null[iter]) {
          name_vec <- c(name_vec, names(zoolist)[iter])
        } else {
          name_vec <- c(name_vec, colnames(zoolist[[iter]]))
        }
      }
      colnames(Y) <- name_vec
    }

    if (all(dim_null)) {
      zoolistfreq <- sapply(Y_in, frequency)
    } else if (all(!dim_null)) {
      zoolistfreq <- sapply(Y_in, frequency)
      zoolistn <- sapply(Y_in, NCOL)
      zoolistfreq <- Reduce(c, mapply(function(x, y) rep(x, each = y), zoolistfreq, zoolistn, SIMPLIFY = FALSE))
    } else {
      zoolistfreq <- c()
      for (iter in 1:length(dim_null)) {
        if (dim_null[iter]) {
          zoolistfreq <- c(zoolistfreq, frequency(Y_in[[iter]]))
        } else {
          zoolistfreq <- c(zoolistfreq, rep(frequency(Y_in[[iter]]), each = ncol(Y_in[[iter]])))
        }
      }
    }
    names(zoolistfreq) <- NULL
    if (all(zoolistfreq %in% c(4, 12))) {
      freq <- ifelse(zoolistfreq == 4, "q", "m")
    } else {
      stop("Only monthly and quarterly frequencies are allowed as time series objects.")
    }
  } else {

  }
  return(list(Y, freq))
}
