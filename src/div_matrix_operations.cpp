#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Computes trace ratio between BTB and Omega
//'
//' @title Computes trace ratio between BTB and Omega
//' @name calculate_r
//'
//' @param BTB A \code{matrix}. Must have dimensions \eqn{(n+4) \times (n+4)}, where
//'        \eqn{n} is the total number of unique observations.
//' @param Omega A \code{matrix}. Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' The trace ratio between BTB and Omega, \eqn{r}, is calculated in the following manner:
//' \enumerate{
//'   \item Calculates the trace of \code{BTB} excluding the first \eqn{2} and last \eqn{3} elements. We call this \eqn{T_1}.
//'   \item Calculates the trace of \code{Omega} excluding the first \eqn{2} and last \eqn{3} elements. We call this \eqn{T_2}.
//'   \item Calculate and return \eqn{r = T_1 / T_2}
//' }
//' \eqn{r} is used in the conversion between \code{sp} and \code{lambda}.
//'
//' @return A \code{double}. The trace ratio between BTB and Omega.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
double calculate_r(const arma::mat& BTB, const arma::mat& Omega) {

  // Get diagonal elements efficiently
  arma::vec diag_BTB = BTB.diag();
  arma::vec diag_Omega = Omega.diag();

  // Create index range excluding first 2 and last 2 elements
  arma::uvec idx = arma::regspace<arma::uvec>(2, diag_BTB.n_elem - 4);

  // Calculate sums of selected diagonal elements
  double sum_BTB = arma::sum(diag_BTB.elem(idx));
  double sum_Omega = arma::sum(diag_Omega.elem(idx));

  return sum_BTB / sum_Omega;
}

//' Computes the inverse matrix of \eqn{B^{T}WB + \lambda \Omega}
//'
//' @title Computes the inverse matrix of \eqn{B^{T}WB + \lambda \Omega}
//' @name calculate_Q
//'
//' @param lambda A \code{double}. The computitional smoothing parameter.
//' @param BTB A \code{matrix}. The inner matrix product \eqn{B^{T}WB}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//' @param Omega A \code{matrix}. The penalty matrix.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' Calculates \eqn{Q = (B^{T}WB + \lambda \Omega)^{-1}}.
//'
//' @return A \eqn{(n+4) \times (n+4)} \code{matrix}, which we denote \eqn{Q}.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
arma::mat calculate_Q(double lambda,
                      const arma::mat& BTB,
                      const arma::mat& Omega) {

  // Calculate BTB + lambda * Omega
  arma::mat M = BTB + lambda * Omega;

  // Solve the system using Cholesky decomposition for better numerical stability
  return arma::solve(M, arma::eye(arma::size(M)), arma::solve_opts::fast);
}

//' Calculates the Smoother Matrix \eqn{S}
//'
//' @title Calculates the Smoother Matrix \eqn{S}
//' @name calculate_S
//'
//' @param lambda A \code{double}. The computitional smoothing parameter.
//' @param weights A \code{numeric} vector of weights. Must be of length \eqn{n}.
//' @param B A \code{matrix}. The B-spline basis matrix.
//'        Must have dimensions \eqn{(n) \times (n+4)}.
//' @param BTB A \code{matrix}. The inner matrix product \eqn{B^{T}WB}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//' @param Omega A \code{matrix}. The penalty matrix.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' Calculates the smoother matrix by \eqn{S = B(B^{T}WB + \lambda \Omega)^{-1}B^{T}}.
//'
//' @return A \eqn{(n) \times (n)} \code{matrix}, which we denote \eqn{S}.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
arma::mat calculate_S(double lambda,
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

//' Calculates the Smoother Matrix \eqn{S}
//'
//' @title Calculates the Smoother Matrix \eqn{S}
//' @name calculate_S2
//'
//' @param weights A \code{numeric} vector of weights. Must be of length \eqn{n}.
//' @param B A \code{matrix}. The B-spline basis matrix.
//'        Must have dimensions \eqn{(n) \times (n+4)}.
//' @param Q A \code{matrix}. The matrix \eqn{Q = (B^{T}WB + \lambda \Omega)^{-1}}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' Calculates the smoother matrix by \eqn{S = BQB^{T}}.
//'
//' @return A \eqn{(n) \times (n)} \code{matrix}, which we denote \eqn{S}.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
arma::mat calculate_S2(const arma::vec& weights,
                       const arma::mat& B,
                       const arma::mat& Q) {

  // Return Smoother matrix
  return B * Q * B.t() * arma::diagmat(weights);
}

//' Calculates the matrix part, \eqn{L}, of the covariance matrix \eqn{\hat{\sigma}^2 \cdot L}
//'
//' @title Calculates the matrix part, \eqn{L} of the covariance matrix \eqn{\hat{\sigma}^2 \cdot L}
//' @name calculate_cov_beta
//'
//' @param BTB A \code{matrix}. The inner matrix product \eqn{B^{T}WB}.
//'        Alternatively, The inner matrix product \eqn{B^{T}WWB}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @param Q A \code{matrix}. The matrix \eqn{Q = (B^{T}WB + \lambda \Omega)^{-1}}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' Calculates \eqn{L} in the covariance matrix of the smoothing spline coefficients
//'
//' \eqn{\mathrm{Cov}[\hat{\beta}] = \hat{\sigma}^2 \cdot L}
//'
//' Where \eqn{L = QB^{T}WBQ}.
//'
//' @return A \eqn{(n+4) \times (n+4)} \code{matrix}, which we denote \eqn{L}.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
arma::mat calculate_cov_beta(const arma::mat& BTB,
                             const arma::mat& Q) {

  // Return Smoother matrix un-scaled cov_beta
  // return arma::symmatu(Q * BTB * Q);
  return Q * BTB * Q;
}

//' Calculates the variance of predicted values
//'
//' @title Calculates the variance of predicted values
//' @name calculate_pred_var
//'
//' @param B_new A \code{matrix}. The possibly new B-spline basis matrix or its derivative.
//'        Must have dimensions \eqn{(m) \times (n+4)}.
//' @param cov_beta A \code{matrix}. The matrix \eqn{\hat{\sigma}^2 \cdot L}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}. See details
//'
//' @details
//' Calculates estimated variances of smoothing spline predicted values.
//' In fact we calculate the diagonal elements of
//'
//' \eqn{\mathrm{Cov}[\hat{\beta}] = \hat{\sigma}^2 \cdot L}
//'
//' Where \eqn{L = QB^{T}W_{a}BQ}, where
//' \itemize{
//'   \item \eqn{Q = (B^{T}WB + \lambda \Omega)^{-1}}
//'   \item \eqn{W_{a}} is either \eqn{W} or \eqn{WW^{T}}
//'   \item \eqn{B} is the original \eqn{n \times (n+4)} B-spline matrix.
//'   \item \eqn{\hat{\sigma}^2 = \mathrm{WRSS} / (n - \mathrm{df})}
//' }
//'
//' @return A \code{numeric} vector of length \eqn{m}. The estimated variances of
//'         the predicted values based on new observations.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
arma::vec calculate_pred_var(const arma::mat& B_new,
                             const arma::mat& cov_beta) {

 // Calculate diagonal elements directly without computing full matrix
 return sum((B_new * cov_beta) % B_new, 1);
}

//' Calculates the effective degrees of freedom \eqn{\mathrm{df}}
//'
//' @title Calculates the effective degrees of freedom \eqn{\mathrm{df}}
//' @name calculate_df
//'
//' @param lambda A \code{double}. The computitional smoothing parameter.
//' @param weights A \code{numeric} vector of weights. Must be of length \eqn{n}.
//' @param B A \code{matrix}. The B-spline basis matrix.
//'        Must have dimensions \eqn{(n) \times (n+4)}.
//' @param BTB A \code{matrix}. The inner matrix product \eqn{B^{T}WB}.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//' @param Omega A \code{matrix}. The penalty matrix.
//'        Must have dimensions \eqn{(n+4) \times (n+4)}.
//'
//' @details
//' Calculates the effective degrees of freedom \eqn{\mathrm{df}} by \eqn{\mathrm{tr}[B(B^{T}WB + \lambda \Omega)^{-1}B^{T}]}.
//' In other words, \eqn{\mathrm{df}} is just the trace of the smoother matrix \eqn{S}.
//'
//' @return A \code{double}. The effective degrees of freedom for the smoothing spline.
//' @export
//'
//' @examples print(1)

// [[Rcpp::export]]
double calculate_df(double lambda,
                    const arma::vec& weights,
                    const arma::mat& B,
                    const arma::mat& BTB,
                    const arma::mat& Omega) {

  // Calculate BTB + lambda * Omega
  arma::mat M = BTB + lambda * Omega;

  // Solve the system using Cholesky decomposition for better numerical stability
  arma::mat S = B * arma::solve(M, B.t() * arma::diagmat(weights),
                                arma::solve_opts::fast);

  double df = arma::trace(S);

  // Return Smoother matrix
  return df;
}
