% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calculate_S}
\alias{calculate_S}
\title{Calculates the Smoother Matrix \eqn{S}}
\arguments{
\item{lambda}{A \code{double}. The computitional smoothing parameter.}

\item{weights}{A \code{numeric} vector of weights. Must be of length \eqn{n}.}

\item{B}{A \code{matrix}. The B-spline basis matrix.
Must have dimensions \eqn{(n) \times (n+4)}.}

\item{BTB}{A \code{matrix}. The inner matrix product \eqn{B^{T}WB}.
Must have dimensions \eqn{(n+4) \times (n+4)}.}

\item{Omega}{A \code{matrix}. The penalty matrix.
Must have dimensions \eqn{(n+4) \times (n+4)}.}
}
\value{
A \eqn{(n) \times (n)} \code{matrix}, which we denote \eqn{S}.
}
\description{
Calculates the Smoother Matrix \eqn{S}
}
\details{
Calculates the smoother matrix by \eqn{S = B(B^{T}WB + \lambda \Omega)^{-1}B^{T}}.
}
\examples{
print(1)
}
