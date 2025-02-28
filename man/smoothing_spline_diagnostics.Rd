% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothing_spline_diagnostics.R
\name{smoothing_spline_diagnostics}
\alias{smoothing_spline_diagnostics}
\title{Smoothing Spline Diagnostics}
\usage{
smoothing_spline_diagnostics(
  data,
  only_optimal_dfs = FALSE,
  weighted = FALSE,
  exclude_not_restricted = TRUE,
  na_rm = TRUE
)
}
\arguments{
\item{data}{A dataset}

\item{only_optimal_dfs}{A \code{logical} value. If \code{TRUE}, only the
optimal values of \eqn{\mathrm{df}(c)} are returned.
Not plots and no tables contaning relationships between
\eqn{\mathrm{df}(c)} and \eqn{\hat{\zeta}} and cross-validation
scores.}

\item{weighted}{A \code{logical} value. If \code{TRUE}, include weighted results.}

\item{exclude_not_restricted}{A \code{logical} value. If \code{FALSE}, exclude
non-restricted results.}

\item{na_rm}{A \code{logical} value. If \code{TRUE}, \code{NA} values are
removed before the diagnostics are performed.}
}
\value{
Returns a \code{list}:
\itemize{
 \item \code{item1:} bla bla
 \item \code{item2:} bla bla
 \item \code{item3:} bla bla
}
}
\description{
Perform diagnostics of a smoothing spline fit based.
}
\details{
Before utilizing the \code{smoothing_spline()} function to fit a smoothing spline
to the \code{data}, it may be advantegous to consider some diagnostics of the
regarded dataset. This function return key diagnostic components such as
\itemize{
 \item \eqn{\mathrm{df}(c)} versus Cross-Validation scores
 \item \eqn{\mathrm{df}(c)} versus \eqn{\hat{\zeta}}
 \item Optimal values of \eqn{\mathrm{df}(c)} based on different variants of
       Generalized Cross-Validation (GCV) and Leave-One-Out (LOO) Cross-Validation
}
These diagnostic components may aid in choosing favorable values of \eqn{\mathrm{df}(c)}
based on \code{data}.
}
\examples{
print(1)
}
