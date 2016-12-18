\name{bvar_data}
\docType{data}
\alias{full_tbl}
\alias{mf_list}
\title{Data}
\description{
  Data to use with the MFBVAR
}
\usage{
  full_tbl
  mf_list
  }
\format{\code{full_tbl} is a \code{\link[tibble]{tibble}} with two columns:
\describe{
    \item{fcst_date}{The date at which the forecast is made.}
  \item{data}{A list-column in which the elements are tibbles containing the data. The variabels are untransformed.}
}
\code{mf_list} is essentially the same as \code{full_tbl}, but with the data appropriately transformed. It is also not a tibble, but a list with elements:
\describe{
    \item{fcst_date}{The date at which the forecast is made.}
  \item{data}{A list in which the elements are tibbles containing the data. The variabels are transformed.}
}
}


