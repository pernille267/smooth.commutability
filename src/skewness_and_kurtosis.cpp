#include <Rcpp.h>
using namespace Rcpp;

//' Calculate sample skewness based on a numeric vector \code{x}.
//'
//' @title Calculate Sample Skewness of a Random Sample
//' @name skewness
//'
//' @param x A \code{numeric} vector that is a random sample.
//' @param na_rm A \code{logical} value. If \code{TRUE}, \code{NA}-values are
//'        removed prior to calculation of sample skewness.
//'
//' @details
//' Calculates the sample skewness of a random sample \eqn{\lbrace x_i \rbrace_{i=1}^{N}}.
//' The sample skewness, attempts to estimate the theoretical skewness, \eqn{\gamma}, by using
//' the following estimator
//'
//' \eqn{\hat{\gamma} = \frac{\sqrt{N(N-1)}}{N-2} \cdot \frac{\frac{1}{N}\sum_{i=1}^{N}(x_i - \overline{x})^3}{\Big[\frac{1}{N}\sum_{i=1}^{N}(x_i - \overline{x})^2\Big]^{1.5}}}.
//'
//' Note that this estimator will be biased if not \eqn{x_i \sim \mathrm{N}(\mu, \sigma^2)}.
//'
//' @return A \code{double} that is the calculated sample skewness.
//'
//' @examples \dontrun{
//'   y <- rlnorm(n = 1000, meanlog = 0, sdlog = 0.25)
//'   skew_y <- skewness(y)
//'   print(skew_y)
//' }
// [[Rcpp::export]]
double skewness(NumericVector x, bool na_rm = true) {
  int n = x.size();
  NumericVector x_clean;

  // Handle NA values
  if (na_rm) {
    x_clean = x[!is_na(x)];
  } else {
    // If any NA present and na_rm is false, return NA
    if (any(is_na(x))) {
      return NA_REAL;
    }
    x_clean = x;
  }

  // Get new size after NA removal
  n = x_clean.size();

  // Check if we have enough data points
  if (n < 3) {
    return NA_REAL;
  }

  double mean_x = mean(x_clean);
  double s2 = 0.0, s3 = 0.0;

  for(int i = 0; i < n; ++i) {
    double diff = x_clean[i] - mean_x;
    s2 += diff * diff;
    s3 += diff * diff * diff;
  }

  s2 = s2 / n;
  s3 = s3 / n;

  // Calculate skewness
  double skew = (s3 / pow(s2, 1.5)) * sqrt((double)n * (n - 1)) / (n - 2);

  return skew;
}

//' Calculate sample excess kurtosis based on a numeric vector \code{x}.
//'
//' @title Calculate Sample Excess Kurtosis of a Random Sample
//' @name kurtosis
//'
//' @param x A \code{numeric} vector that is a random sample.
//' @param na_rm A \code{logical} value. If \code{TRUE}, \code{NA}-values are
//'        removed prior to calculation of sample excess kurtosis.
//'
//' @details
//' Calculates the sample excess kurtosis of a random sample \eqn{\lbrace x_i \rbrace_{i=1}^{N}}.
//' The sample excess kurtosis, attempts to estimate the theoretical excess kurtosis, \eqn{\Kappa}, by using
//' the following estimator
//'
//' \eqn{\hat{\Kappa} = \frac{N(N+1)}{(N-1)(N-2)(N-3)} \cdot \frac{\sum_{i=1}^{N}(x_i - \overline{x})^4}{s^4} - 3 \cdot \frac{(N-1)^2}{(N-2)(N-3)}}.
//'
//' Here, \eqn{s^2} is the unbiased sample variance.
//' Note that this estimator will be biased if not \eqn{x_i \sim \mathrm{N}(\mu, \sigma^2)}.
//' Note also that \eqn{\hat{\Kappa}} estimates the theoretical excess kurtosis and not the raw kurtosis.
//' Therefore, one should add \eqn{3} to the output to get the estimated raw kurtosis.
//'
//' @return A \code{double} that is the calculated sample excess kurtosis.
//'
//' @examples \dontrun{
//'   y <- rnorm(n = 1000)
//'   kurt_y <- kurtosis(y)
//'   print(kurt_y)
//' }
// [[Rcpp::export]]
double kurtosis(NumericVector x, bool na_rm = true) {
  int n = x.size();
  NumericVector x_clean;

  // Handle NA values
  if (na_rm) {
    x_clean = x[!is_na(x)];
  } else {
    // If any NA present and na_rm is false, return NA
    if (any(is_na(x))) {
      return NA_REAL;
    }
    x_clean = x;
  }

  // Get new size after NA removal
  n = x_clean.size();

  // Check if we have enough data points
  if (n < 3) {
    return NA_REAL;
  }

  double mean_x = mean(x_clean);
  double s2 = 0.0, s4 = 0.0;

  for(int i = 0; i < n; ++i) {
    double diff = x_clean[i] - mean_x;
    s2 += diff * diff;
    s4 += diff * diff * diff * diff;
  }

  s2 = s2 / (n - 1);

  double fact1 = (double)n * (n + 1) / (n - 1) / (n - 2) / (n - 3);
  double fact2 = 3.0 * pow(n - 1, 2) / (n - 2) / (n - 3);

  // Calculate excess kurtosis
  double kurt = fact1 * (s4 / pow(s2, 2.0)) - fact2;

  return kurt;
}

