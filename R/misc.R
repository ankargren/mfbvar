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

all_equal <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}


