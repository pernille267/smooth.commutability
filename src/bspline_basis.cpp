#include <Rcpp.h>
using namespace Rcpp;

// Derive the B-spline basis function value at a value x
//'
//' @title Derive the B-spline basis function value at a value x
//' @name bspline_basis
//'
//' @param j A \code{integer} referring to the knot interval where \code{x} is to be evaluated.
//' @param m A \code{integer} signifying the degree of the B-spline basis.
//' @param knots A \code{numeric} vector containing the knots for the B-spline basis.
//' @param x A \code{double} corresponding to the value which we seek the B-spline function value.
//'
//' @description Computes the B-spline function value evaluated at \code{x} for the \code{j}-th knot, when the degree of the B-spline basis is \code{m}.
//'
//' @details Used to recursively build the B-spline basis matrix B.
//'
//' @return A \code{double} that is the function value of the B-spline function evaluated at \code{x} for the \code{j}-th knot given the B-spline basis degree \code{m}.
//'
//' @examples \dontrun{
//'   print(1)
//' }


// [[Rcpp::export]]
double bspline_basis(int j, int m, NumericVector knots, double x) {
  // When degree m = 0, return 0 if x is in knot interval
  if(m == 0){
    if(knots[j] == knots[j + 1]){
      return 0.0;
    }
    else if(knots[j] <= x && x < knots[j + 1]){
      return 1.0;
    }
    else{
      return 0.0;
    }
  }
  // When degree m > 0, do the recursion
  else{
    // First term
    double basis_left = 0.0;

    // Second term
    double basis_right = 0.0;

    // Evaluate the first term
    if(knots[j + m] != knots[j]){
      basis_left = (x - knots[j]) / (knots[j + m] - knots[j]) * bspline_basis(j, m - 1, knots, x);
    }

    // Evaluate the second term
    if(knots[j + m + 1] != knots[j + 1]){
      basis_right = (knots[j + m + 1] - x) / (knots[j + m + 1] - knots[j + 1]) * bspline_basis(j + 1, m - 1, knots, x);
    }

    return basis_left + basis_right;
  }
}

// Derive the B-spline basis matrix given values of x and the degree m
//'
//' @title Derive the B-spline basis matrix given values of x and the degree m
//' @name bspline_basis_matrix
//'
//' @param x A \code{numeric} vector containing the predictior values.
//' @param knots A \code{numeric} vector containing the knots for the B-spline basis.
//' @param m An \code{integer} that signify the desired degree of the B-spline basis. Default is 3 (cubic).
//'
//' @description Computes the B-spline basis matrix given the \code{knots} and the desired degree, \code{m}.
//'
//' @details Recursively build the B-spline basis matrix B. The construction of B uses the Cox-de Boor recursion formula.
//'
//' @return The constructed B-spline basis \code{matrix}.
//'
//' @examples \dontrun{
//'   print(1)
//' }


// [[Rcpp::export]]
NumericMatrix bspline_basis_matrix(NumericVector x, NumericVector knots, int m = 3){

  // Handling boundary knots
  NumericVector boundary_knots(2);
  boundary_knots[0] = min(knots);
  boundary_knots[1] = max(knots);
  NumericVector additional_knots(2 * (m + 1));

  // Create M equal boundary knots on both left and right side of the interval
  int ind = 0;
  for(int i = 0; i < additional_knots.size() / 2; ++i){
    additional_knots[ind] = boundary_knots[0];
    additional_knots[ind + 1] = boundary_knots[1];
    ind+=2;
  }

  // Gather interior knots and replications of boundary knots
  NumericVector augmented_knots(knots.size() + 2 * (m + 1));
  for(int i = 0; i < 2 * (m + 1); ++i){
    augmented_knots[i] = additional_knots[i];
  }
  for(int i = ind; i < augmented_knots.size(); ++i){
    augmented_knots[i] = knots[i - ind];
  }
  NumericVector sorted_augmented_knots = clone(augmented_knots);
  NumericVector sorted_x = clone(x);

  std::sort(sorted_x.begin(), sorted_x.end());
  std::sort(sorted_augmented_knots.begin(), sorted_augmented_knots.end());

  int n_intervals = sorted_augmented_knots.size() - m - 1;
  int n_points = x.size();

  NumericMatrix basis_matrix(n_points, n_intervals);

  for (int i = 0; i < n_points; i++){
    for (int j = 0; j < n_intervals; j++){
      basis_matrix(i, j) = bspline_basis(j, m, sorted_augmented_knots, sorted_x[i]);
    }
  }

  return basis_matrix;
}
