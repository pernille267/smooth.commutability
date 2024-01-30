#' Convert the desired effective degrees of freedom to the computational smoothing parameter lambda
#'
#' @param df A \code{double} that corresponds to the desired effective degrees of freedom for the smoothing spline model.
#' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of the predictor values, \code{x}.
#' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
#' @param Omega The O'Sullivan penalty \code{matrix} for the B-spline basis.
#'
#' @return A \code{double} that is the associated value of lambda to the provided \code{df}.
#' @export
#'
#' @examples print(1)

to_lambda <- function(df, B, BTB, Omega){
  # Numerical optimization of the absolute difference between Trace(S) and the target df
  optimal_error_df <- optimise(f = error_df, lower = -1.5, upper = 1.5, df = df, B = B, BTB = BTB, Omega = Omega)
  lambda <- optimal_error_df$minimum
  # Convert scale-free smoothing parameter (aka "spar") to computational lambda
  r <- sum(diag(BTB)) / sum(diag(Omega))
  clambda <- r * 256**(3 * lambda - 1)
  return(clambda)
}

