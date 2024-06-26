% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothing_spline.R
\name{smoothing_spline}
\alias{smoothing_spline}
\title{Fit a Smoothing Spline Model with O'Sullivan Penalty}
\usage{
smoothing_spline(data, df = NULL, lambda = NULL, df_max = 7.5)
}
\arguments{
\item{data}{A \code{data.table}, \code{data.frame} or \code{list} containing the measurements.
It should include two specific variables named 'MP_B' and 'MP_A'.}

\item{df}{Optional; a \code{numeric} value specifying the degrees of freedom for the spline model.
Should be between 2 and the number of clinical samples in \code{data}.
The default value (\code{NULL}) triggers automatic selection using n-fold cross-validation.}

\item{lambda}{Optional; a non-negative \code{numeric} value specifying the penalty strength.
The default (\code{NULL}) results in the function choosing an optimal value based on
n-fold cross-validation to minimize the mean squared prediction error.}

\item{df_max}{A \code{numeric} value between 2 and the number of clinical samples in \code{data}.
If n-fold cross-validation is used to obtain the optimal effective number of degrees of freedom, what is its upper limit.}
}
\value{
A \code{list} containing the components of the smoothing spline fit, including:
        - \code{fitted}: Fitted values of the spline model.
        - \code{residuals}: Residuals from the fitted model.
        - \code{df}: Effective degrees of freedom used.
        - \code{lambda}: Penalty parameter used.
        - \code{cv_crit}: Criterion value from cross-validation.
        - Additional model components.
        These components are often utilized in subsequent analyses and visualizations.
}
\description{
This function fits a smoothing spline model to the provided data using an O'Sullivan penalty matrix. It is designed to work with various data types and offers flexibility in model tuning through degrees of freedom (`df`) and penalty parameter (`lambda`).
}
\examples{
# Example using a dataset with columns MP_B and MP_A
library(data.table)
data_example <- data.table(MP_B = 1:10, MP_A = rnorm(10))
result <- smoothing_spline(data_example, df = 5)
print(result$fitted)
}
