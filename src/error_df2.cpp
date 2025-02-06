#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Computes the absolute effective degrees of freedom error
//'
//' @title Computes the absolute effective degrees of freedom error
//' @name error_df2
//'
//' @param sp A \code{double}, and is a scale-free variant of the computational lambda. Given this value of \code{lambda}, the trace of the smoother matrix is calculated and the absolute deviation between this and the given \code{df} is recorded.
//' @param df A \code{double} that corresponds to the target effective degrees of freedom for the smoothing spline model.
//' @param weights A \code{numeric} vector sorted according to the order of the predictor \code{x}.
//' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of the predictor values, \code{x}.
//' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
//' @param Omega The O'Sullivan penalty \code{matrix} for the B-spline basis.
//' @param r The trace ratio of the matrices
//'
//' @description
//' Given a value of \code{lambda}, the corresponding effective degrees of freedom is calculated. This value of the degrees of freedom is compared to \code{df}. This comparison is done by taking the absolute difference between the effective degrees of freedom associated with the \code{lambda} given and the provided value of \code{df}. This function is primarily used to convert a value of \code{df} to a suitable value of \code{lambda} using numerical optimization. See \code{to_lambda} for more details.
//'
//'
//' @return A \code{double} that signify the absolute deviation between the trace of the smoother matrix (given the value of \code{lambda}) and the provided value of \code{df}.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
double error_df2(double sp,
                 double df,
                 const arma::vec& weights,
                 const arma::mat& B,
                 const arma::mat& BTB,
                 const arma::mat& Omega,
                 double r = 1.0) {
  // Convert scale-free smoothing parameter to computational lambda
  double lambda = r * std::pow(256.0, 3.0 * sp - 1.0);

  // Calculate BTB + lambda * Omega
  arma::mat M = BTB + lambda * Omega;

  // Solve the system using Cholesky decomposition for better numerical stability
  arma::mat S = B * arma::solve(M, B.t() * arma::diagmat(weights),
                                arma::solve_opts::fast);

  // Calculate trace efficiently using diagonal elements
  double candidate_df = arma::trace(S);

  // Return absolute difference
  return std::abs(candidate_df - df);
}
