#' Convert the desired effective degrees of freedom to the computational smoothing parameter lambda
#'
#' @param df A \code{double} that corresponds to the desired effective degrees of freedom for the smoothing spline model.
#' @param weights A \code{numeric} vector sorted according to the order of the predictor \code{x}.
#' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of the predictor values, \code{x}.
#' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
#' @param Omega The O'Sullivan penalty \code{matrix} for the B-spline basis.
#' @param r The trace ratio of the matrices
#'
#' @return A \code{double} that is the associated value of lambda to the provided \code{df}.
#' @export
#'
#' @examples print(1)

to_lambda <- function(df, weights, B, BTB, Omega, r = 1){
  optimal_error_df <- optimise(f = error_df, lower = -1.5, upper = 1.5, df = df, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
  if(optimal_error_df$objective > 0.05){
    optimal_error_df <- optimise(f = error_df, lower = -2.0, upper = 3.0, df = df, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
    if(optimal_error_df$objective > 0.05){
      return(NA)
    }
  }
  sp <- optimal_error_df$minimum
  lambda <- r * 256**(3 * sp - 1)
  return(lambda)
}

