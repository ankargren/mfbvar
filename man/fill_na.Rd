\name{fill_na}
\alias{fill_na}
\title{Fills NAs with the next non-NA value}
\description{
The function fills elements with \code{NA} with the next non-\code{NA} value (so that quarterly averages observed at the end of the quarter are assumed as observations for the remaining months of the quarter).
}
\usage{
fill_na(Y)
}
\arguments{
\item{Y}{The data matrix of size \code{n_T * n_vars} containing some missing values (\code{NA}).}
}
\value{
\item{x}{A matrix with no \code{NA}s.}
}

