% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fun_of_replicates2}
\alias{fun_of_replicates2}
\title{Calculate a Statistic for Each Sample}
\usage{
fun_of_replicates2(data, fun = "mean")
}
\arguments{
\item{data}{A \code{data.table} or \code{list} object. Must contain
\code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.}

\item{fun}{A \code{character} string. Which statistic is to be calculated for
each SampleID. Possible choices include \code{mean}, \code{var} (variance),
\code{sd} (standard deviaton), \code{cv} (coefficient of variation),
\code{median}, \code{min} (minimum) and \code{max} (maximum).}
}
\value{
A \code{list} with elements \code{SampleID}, \code{MP_A} and \code{MP_B}.
        \code{MP_A} and \code{MP_B} contains the sample-wise statistics.
}
\description{
Calculates a chosen statistic over the replicates.
}
\details{
Calculate a Statistic for Each Sample

This function handles NA-values automatically.
}
\examples{
\dontrun{
  print(1)
}
}
