#include <Rcpp.h>
using namespace Rcpp;

// Obtain optimal upper limit for constrained search grid for optimal df
//'
//' @title Obtain optimal upper limit for constrained search grid for optimal df
//' @name obtain_df_max
//'
//' @param df A \code{numeric} vector signifying the grid in which second derivatives of zeta are evaluated
//' @param m A \code{numeric} vector signifying the evaluated second derivatives of zeta
//' @param threshold A \code{double}. Which relative value of the maximum found second derivative should be used as tolerance. Must be between 0 and 1.
//'
//' @description Computes the first value of \code{df} where the curve of zeta starts to increase after initial decrease.
//'
//' @details Used to constrain the search grid for an optimal value of \code{df}.
//'
//' @return A \code{double}. The upper limit for the constrained search grid for optimal \code{df}.
//'
//' @examples \dontrun{
//'   print(1)
//' }

// [[Rcpp::export]]
double obtain_df_max(NumericVector df, NumericVector second_deriv, double threshold = 0.05) {
  int max_second_deriv = which_max(second_deriv);
  if (df.length() != second_deriv.length()) {
    stop("df and second_deriv must have the same length");
  }

  for (int i = 0; i < df.length(); i++) {
    if (df[i] > df[max_second_deriv] && second_deriv[i] < threshold * second_deriv[max_second_deriv]) {
      return df[i];
    }
  }

  // If no value is found, return NA
  return df[max_second_deriv];
}
