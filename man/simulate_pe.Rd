% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_mspe.R
\name{simulate_pe}
\alias{simulate_pe}
\title{Simulate Smoothing Spline Prediction Errors}
\usage{
simulate_pe(
  parameters,
  m = 1,
  extrapolate_ok = FALSE,
  shift = FALSE,
  attempt_fast = FALSE
)
}
\arguments{
\item{parameters}{A \code{list}, \code{data.table} or \code{data.frame}.
The parameters and corresponding values to be used in
the simulation. Either \code{df} or \code{df_max}, must
be included. See \code{sim_eqa_data()} for other
parameters.}

\item{m}{An \code{integer}. Must be larger than or equal to \code{1L}.}

\item{extrapolate_ok}{A non-missing \code{logical} value. If \code{TRUE},
simulated evaluated material data are allowed to
produce measurements beyond the support of the
clinical samples IVD-MD measurements. Defaults to
\code{FALSE}.}

\item{shift}{A non-missing \code{logical} value. If \code{TRUE}, the roles
of \code{MP_A} and \code{MP_B} are allowed to shift. Can only
happen if \code{lambda < 0.5}. Defaults to \code{FALSE}.}

\item{attempt_fast}{A non-missing \code{logical} value. If \code{TRUE}, the
\code{smooth.spline()} is used to make predictions.}
}
\value{
A \code{numeric} vector of length \code{m}. The simulated prediction errors.
}
\description{
Simulates smoothing spline prediction errors for each of the \code{m}
evaluated materials.
}
\details{
Simulate a clinical sample dataset and an evaluated material dataset. The
evaluated material are assumed to be commutable and will thus belong to the
same population as the clinical samples. Then the corresponding
prediction error(s) is(are) calculated. The prediction errors is of course
calculated by the arithmetic difference between the observed value(s) of
the evaluated material and the predicted value using the smoothing spline
model.

Note: This function is generally not used directly, but as part of the
more comprehensive function \code{simulate_mspe()}.

This function will never throw an error. If the invalid arguments are passed
the output will instead produce \code{NA}.

General note: Any function named \code{simulate_***}, is generally not
intended to be used by end-users. They are used in research and validity
testing.
}
\examples{
# Required packages
library(fasteqa)

# Used parameters
parameters <- list(n = 25,
                   R = 3,
                   cvx = 0.01,
                   cvy = 0.01,
                   obs_tau = 6,
                   df = 6,
                   type = 2)

# Simulate 1000 smoothing spline prediction errors at y = 6
ss_pes <- replicate(n = 1000,
                    expr = simulate_pe(parameters),
                    simplify = TRUE)

# Calculate prediction bias at y = 6
print(round(mean(ss_pes, na.rm = TRUE), 4L))

# Calculate variance of prediction error at y = 6
print(round(var(ss_pes, na.rm = TRUE), 4L))


}
