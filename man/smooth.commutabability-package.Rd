\name{smooth.commutabability-package}
\alias{smooth.commutabability}
\docType{package}
\title{
  Commutability Analysis for Clinical Laboratory Data Using Smoothing Splines
}
\description{
  Provides methods for analyzing and visualizing commutability of clinical
  data using spline models. Implements
  novel algorithms for commutability evaluation.
}
\details{
  Key functions:
  \itemize{
    \item \code{\link{smoothing_spline}()}: Fitting a smoothing spline
    \item \code{\link{plot.smoothing_spline}()}: Visualization of model fits
    \item \code{\link{predict_smoothing_spline}()}: Construct smoothing spline prediction intervals
    \item \code{\link{estimate_zeta_ss}()}: Estimate Differences in nonselectivity
  }
}
\author{
Pernille Kjeilen Fauskanger <pernille.fauskanger@noklus.no>
}
\references{
  Wahba, G. (1983). Bayesian "Confidence Intervals" for the Cross-Validated Smoothing Spline. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology} 45(1): 133-150.
  \url{https://doi.org/10.1111/j.2517-6161.1983.tb01239.x}

  Silverman, B. W. (1985). Some Aspects of the Spline Smoothing Approach to Non-Parametric Regression Curve Fitting. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 47(1), 1-21.
  \url{https://doi.org/10.1111/j.2517-6161.1985.tb01327.x}

  Kim, Y.-J., and Gu, C. (2004). Smoothing Spline Gaussian Regression: More Scalable Computation via Efficient Approximation. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, 66(2): 337-356.
  \url{https://doi.org/10.1046/j.1369-7412.2003.05316.x}

  Hastie, T., Tibshirani, R., & Friedman, J. (2009). The Elements of Statistical Learning. \emph{Springer New York}.
  \url{https://doi.org/10.1007/978-0-387-84858-7}

  Wand, M.P. and Ormerod, J.T. (2008), ON SEMIPARAMETRIC REGRESSION WITH O'SULLIVAN PENALIZED SPLINES \emph{Australian & New Zealand Journal of Statistics} 50: 179-198.
  \url{https://doi.org/10.1111/j.1467-842X.2008.00507.x}

  Sandberg S. et al. (2023) for the IFCC Working Group on Commutability in Metrological Traceability, Recommendations for Setting a Criterion and Assessing Commutability of Sample Materials Used in External Quality Assessment/Proficiency Testing Schemes, \emph{Clinical Chemistry}. 69(11): 1227-1237.
  \url{https://doi.org/10.1093/clinchem/hvad135}

  Fauskanger P.K., et al. (2025) Quantification of Difference in Nonselectivity Between In Vitro Diagnostic Medical Devices. \emph{Biometrical Journal}. 67: e70032.
  \url{https://doi.org/10.1002/bimj.70032}
}
\keyword{ package }
\seealso{
  Useful links:
  \itemize{
    \item \code{\link{smoothing_spline}()}
  }
}
\examples{
  # Basic commutability analysis
  print(1)
}
