% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{resample_imprecision}
\alias{resample_imprecision}
\title{Resample imprecision estimates based on clustered EQA clinical sample data}
\arguments{
\item{data}{A \code{list} or a \code{data.table}. Must contain \code{SampleID},
\code{ReplicateID}, \code{MP_A} and \code{MP_B}. The ID variables \code{SampleID}
and \code{ReplicateID} must be of character type for the function to operate correctly.}
}
\value{
A \code{list} containing the resampled imprecision.
}
\description{
Resample imprecision estimates based on clustered EQA clinical sample data
}
\details{
This function is a very efficient method to resample repeatability measure estimates
based on clinical sample data in \code{data}.
}
\examples{
library(data.table)
fictive_data <- simulate_eqa_data2(list(n = 25, R = 3, cvx = 0.01, cvy = 0.01), AR = TRUE)
impr <- replicate(n = 5, expr = resample_imprecision(fictive_data), simplify = FALSE)
impr <- rbindlist(impr)
print(impr)
}
