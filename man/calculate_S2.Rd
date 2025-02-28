% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calculate_S2}
\alias{calculate_S2}
\title{Calculates the Smoother Matrix \eqn{S}}
\arguments{
\item{weights}{A \code{numeric} vector of weights. Must be of length \eqn{n}.}

\item{B}{A \code{matrix}. The B-spline basis matrix.
Must have dimensions \eqn{(n) \times (n+4)}.}

\item{Q}{A \code{matrix}. The matrix \eqn{Q = (B^{T}WB + \lambda \Omega)^{-1}}.
Must have dimensions \eqn{(n+4) \times (n+4)}.}
}
\value{
A \eqn{(n) \times (n)} \code{matrix}, which we denote \eqn{S}.
}
\description{
Calculates the Smoother Matrix \eqn{S}
}
\details{
Calculates the smoother matrix by \eqn{S = BQB^{T}}.
}
\examples{
print(1)
}
