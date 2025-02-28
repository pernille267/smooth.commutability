% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{obtain_df_max}
\alias{obtain_df_max}
\title{Obtain optimal upper limit for constrained search grid for optimal df}
\arguments{
\item{df}{A \code{numeric} vector signifying the grid in which second derivatives of zeta are evaluated}

\item{m}{A \code{numeric} vector signifying the evaluated second derivatives of zeta}

\item{threshold}{A \code{double}. Which relative value of the maximum found second derivative should be used as tolerance. Must be between 0 and 1.}
}
\value{
A \code{double}. The upper limit for the constrained search grid for optimal \code{df}.
}
\description{
Computes the first value of \code{df} where the curve of zeta starts to increase after initial decrease.
}
\details{
Used to constrain the search grid for an optimal value of \code{df}.
}
\examples{
\dontrun{
  print(1)
}
}
