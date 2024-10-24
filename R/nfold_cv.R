#' N-fold cross-validation estimate of the mean squared prediction error of the smoothing spline model
#'
#' @param sp A non-negative \code{double}. This is a scaled version of the raw smoothing parameter lambda and is usually what is optimized.
#' @param y A \code{numeric} vector sorted according to the order of the predictor \code{x}.
#' @param weights A \code{numeric} vector sorted according to the order of the predictor \code{x}.
#' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of \code{x} .
#' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
#' @param Omega The penalty \code{matrix} for the B-spline basis.
#' @param cross_validation A \code{logical} value. If set to \code{TRUE}, you are effectively prompting that the lambda inputs are equivalent to the spar smoothing parameter.
#' @param r The trace ratio of the matrices
#'
#' @return A \code{double} that signify the n-fold cross-validation score
#' @export
#'
#' @examples print(1)

nfold_cv <- function(sp, y, weights, B, BTB, Omega, cross_validation = FALSE, r = 1){
  if(cross_validation){
    lambda <- r * 256**(3 * sp - 1)
  }
  else{
    lambda <- sp
  }
  S <- B %*% solve(BTB + lambda * Omega) %*% t(B) %*% diag(weights)
  fitted <- S %*% y
  return(mean(weights * ((y - fitted) / (1 - diag(S)))**2))
}

