#include <Rcpp.h>
using namespace Rcpp;

//' Estimate weighted local average curve.
//'
//' @title Estimate Weighted Local-Average Curve
//' @name local_average
//'
//' @param x A \code{numeric} vector. Must contain the stochastic process values
//'        that are used in estimating the weighted local average curve.
//' @param weights A \code{numeric} vector. Must contain weights for each of the stochastic
//'        process values. Set all weights to \code{1} to bypass weighting.
//' @param window A \code{integer} value. The window-width for estimating averages
//'        at each point. See details.
//'
//' @details
//' Calculates the estimated weighted local-averages for the sample \eqn{\lbrace x_i \rbrace_{i=1}^{N}}
//' and weights \eqn{\lbrace w_i \rbrace_{i=1}^{N}}.
//' We assume the following model:
//'
//' \eqn{x_i = f(x_i) + \epsilon_i}.
//'
//' We attempt to estimate the signal \eqn{f(x_i)}, using a weighted local-average model:
//'
//' \eqn{\hat{f}(x_i) = \frac{1}{w_i \cdot (\min\lbrace N, i + \Delta \rbrace - \max\lbrace 1, i - \Delta \rbrace)} \cdot \sum_{j = \max\lbrace 1, i - \Delta \rbrace}^{\min\lbrace N, i + \Delta \rbrace} x_j}
//'
//' Here, \eqn{\Delta} is half of the value given in \code{window}. \code{NA} values are allowed,
//' and will be silently removed if present.
//'
//' @return A \code{numeric} vector of local-average estimates.
//'
//' @examples \dontrun{
//'   print(1)
//' }
// [[Rcpp::export]]
NumericVector local_average(NumericVector x, NumericVector weights, int window) {
  int n = x.length();
  int half_window = window / 2;
  NumericVector result(n);

  // Handle edge cases and middle cases
  for(int i = 0; i < n; i++) {
    double sum = 0.0;
    int count = 0;

    // Define window boundaries
    int start = std::max(0, i - half_window);
    int end = std::min(n - 1, i + half_window);

    // Calculate average for current window
    for(int j = start; j <= end; j++) {
      if (!R_IsNA(x[j]) && !R_IsNA(weights[j])) {
        sum += x[j] * weights[j];
        count++;
      }
    }

    // Store result
    result[i] = (count > 0) ? sum/count : NA_REAL;
  }

  return result;
}
