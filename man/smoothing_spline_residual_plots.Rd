% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/suggest_df.R
\name{smoothing_spline_residual_plots}
\alias{smoothing_spline_residual_plots}
\title{Plot Smoothing Spline Residual Based plots}
\usage{
smoothing_spline_residual_plots(
  smoothing_spline_object,
  which_plot = c("all", "histogram_studentized", "scatter_studentized", "histogram_raw",
    "scatter_raw")
)
}
\arguments{
\item{smoothing_spline_object}{A \code{object} of either class \code{smooth.spline} or \code{smoothing_spline}.}

\item{which_plot}{A \code{character} string determining which type of residual plot to be plotted. Possible entries include \code{all} (for all four plot variants),
\code{histogram_studentized}, \code{scatter_studentized}, \code{histogram_raw} and \code{scatter_raw}.}
}
\value{
A \code{list} of \code{ggplot2} objects if \code{which_plot = "all"}. Otherwise a single \code{ggplot2} object, containing necessary information on the plot.
}
\description{
This function is used to plot residual plots and residual histograms, meant for model diagnosis and validation purposes.
}
\details{
Note that choosing either \code{scatter_studentized} or \code{histogram_studentized} will use studentized residuals where the
         where residuals are scaled by (sig^2 * (1-s_ii))^(-1/2), where s_ii is the leverage for the i-th smallest predictor value.
}
\examples{
# Example using a fictive dataset
sim_pars <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10, cve = 0)
sim_data <- simulate_eqa_data2(sim_pars, type = 2, AR = TRUE)
smoothing_spline_fit <- smooth.spline(x = sim_data$MP_B, y = sim_data$MP_A, df = 5)
smoothing_spline_residual_plots(smoothing_spline_fit, "histogram_studentized")

}
