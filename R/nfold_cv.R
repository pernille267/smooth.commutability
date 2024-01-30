#' N-fold cross-validation estimate of the mean squared prediction error of the smoothing spline model
#'
#' @param lambda A non-negative \code{double}. This value is usually what is optimized.
#' @param y A \code{numeric} vector sorted according to the order of the predictor \code{x}.
#' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of \code{x} .
#' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
#' @param Omega The penalty \code{matrix} for the B-spline basis.
#' @param cross_validation A \code{logical} value. If set to \code{TRUE}, you are effectively prompting that the lambda inputs are equivalent to the spar smoothing parameter.
#'
#' @return A \code{double} that signify the n-fold cross-validation score
#' @export
#'
#' @examples print(1)

nfold_cv <- function(lambda, y, B, BTB, Omega, cross_validation = FALSE){
  if(cross_validation){
    r <- sum(diag(BTB)) / sum(diag(Omega))
    clambda <- r * 256**(3 * lambda - 1)
  }
  else{
    clambda <- lambda
  }
  S <- B %*% solve(BTB + clambda * Omega) %*% t(B)
  fitted <- S %*% y
  return(mean(((y - fitted) / (1 - diag(S)))**2))
}
