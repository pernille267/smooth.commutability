#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Helper function
arma::mat calculate_S1(double lambda,
                       const arma::vec& weights,
                       const arma::mat& B,
                       const arma::mat& BTB,
                       const arma::mat& Omega) {

  // Calculate BTB + lambda * Omega
  arma::mat M = BTB + lambda * Omega;

  // Solve the system using Cholesky decomposition for better numerical stability
  arma::mat S = B * arma::solve(M, B.t() * arma::diagmat(weights),
                                arma::solve_opts::fast);

  // Return Smoother matrix
  return S;
}

//' N-fold cross-validation estimate of the mean squared prediction error of the smoothing spline model
//'
//' @title N-fold cross-validation estimate of the mean squared prediction error of the smoothing spline model
//' @name nfold_cv2
//'
//' @param sp A non-negative \code{double}. This is a scaled version of the raw smoothing parameter lambda and is usually what is optimized.
//' @param y A \code{numeric} vector sorted according to the order of the predictor \code{x}.
//' @param weights A \code{numeric} vector sorted according to the order of the predictor \code{x}.
//' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of \code{x} .
//' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
//' @param Omega The penalty \code{matrix} for the B-spline basis.
//' @param cross_validation A \code{logical} value. If set to \code{TRUE}, you are effectively prompting that the lambda inputs are equivalent to the spar smoothing parameter.
//' @param r The trace ratio of the matrices
//'
//' @return A \code{double} that signify the n-fold cross-validation score
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
double nfold_cv2(double sp,
                 const arma::vec& y,
                 const arma::vec& weights,
                 const arma::mat& B,
                 const arma::mat& BTB,
                 const arma::mat& Omega,
                 bool cross_validation = false,
                 double r = 1.0) {

  double lambda;
  if (cross_validation) {
    lambda = r * std::pow(256.0, 3.0 * sp - 1.0);
  } else {
    lambda = sp;
  }

  arma::mat S = calculate_S1(lambda, weights, B, BTB, Omega);

  // Compute fitted values
  arma::vec fitted = S * y;

  // Extract diagonal elements of S
  arma::vec S_diag = S.diag();

  // Compute CV score
  arma::vec residuals = y - fitted;
  arma::vec cv_terms = weights % arma::square(residuals / (1.0 - S_diag));

  return arma::mean(cv_terms);
}

//' Generalized cross-validation estimate of the mean squared prediction error of the smoothing spline model
//'
//' @title Generalized cross-validation estimate of the expected mean squared prediction error of the smoothing spline model
//' @name gcv2
//'
//' @param sp A non-negative \code{double}. This is a scaled version of the raw smoothing parameter lambda and is usually what is optimized.
//' @param y A \code{numeric} vector sorted according to the order of the predictor \code{x}.
//' @param weights A \code{numeric} vector sorted according to the order of the predictor \code{x}.
//' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of \code{x} .
//' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
//' @param Omega The penalty \code{matrix} for the B-spline basis.
//' @param cross_validation A \code{logical} value. If set to \code{TRUE}, you are effectively prompting that the lambda inputs are equivalent to the spar smoothing parameter.
//' @param r The trace ratio of the matrices
//' @param smudge A \code{double} larger than or equal to 1. Set to \code{1.2} or
//'        \code{1.4} to avoid undersmoothing.
//'
//' @return A \code{double} that signify the GCV cross-validation score.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
double gcv2(double sp,
            const arma::vec& y,
            const arma::vec& weights,
            const arma::mat& B,
            const arma::mat& BTB,
            const arma::mat& Omega,
            bool cross_validation = false,
            double r = 1.0,
            double smudge = 1.4) {

 // Initialize computational smoothing parameter
 double lambda;
  if (cross_validation) {
    lambda = r * std::pow(256.0, 3.0 * sp - 1.0);
  } else {
    lambda = sp;
  }

  // Calculate smoother matrix S
  arma::mat S = calculate_S1(lambda, weights, B, BTB, Omega);

  // Compute fitted values of smoothing spline
  arma::vec fitted_values = S * y;

  // Calculate mean of diagonal elements of S
  double S_diag_mean = arma::mean(smudge * S.diag());

  // Compute residuals
  arma::vec residuals = y - fitted_values;

  // Calculate GCV score
  arma::vec cv_terms = weights % arma::square(residuals / (1.0 - S_diag_mean));

  return arma::mean(cv_terms);
}
