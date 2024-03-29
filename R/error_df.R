#' Computes the absolute effective degrees of freedom error
#'
#' @param lambda A \code{double}, and is a scale-free variant of the computational lambda. Given this value of \code{lambda}, the trace of the smoother matrix is calculated and the absolute deviation between this and the given \code{df} is recorded.
#' @param df A \code{double} that corresponds to the target effective degrees of freedom for the smoothing spline model.
#' @param B A B-spline basis \code{matrix} evaluated at the sorted and unique values of the predictor values, \code{x}.
#' @param BTB A \code{matrix} that is the matrix product between the transpose of \code{B} and \code{B}.
#' @param Omega The O'Sullivan penalty \code{matrix} for the B-spline basis.
#'
#' @description
#' Given a value of \code{lambda}, the corresponding effective degrees of freedom is calculated. This value of the degrees of freedom is compared to \code{df}. This comparison is done by taking the absolute difference between the effective degrees of freedom associated with the \code{lambda} given and the provided value of \code{df}. This function is primarily used to convert a value of \code{df} to a suitable value of \code{lambda} using numerical optimization. See \code{to_lambda} for more details.
#'
#'
#' @return A \code{double} that signify the absolute deviation between the trace of the smoother matrix (given the value of \code{lambda}) and the provided value of \code{df}.
#' @export
#'
#' @examples print(1)

error_df <- function(lambda, df, B, BTB, Omega){
  # Convert scale-free smoothing parameter (aka "spar") to computational lambda
  r <- sum(diag(BTB)) / sum(diag(Omega))
  clambda <- r * 256**(3 * lambda - 1)
  # Calculate the smoother matrix
  S <- B %*% solve(BTB + clambda * Omega) %*% t(B)
  # Calculate the absolute difference between Trace(S) and df
  candidate_df <- sum(diag(S))
  df_error <- abs(candidate_df - df)
  return(df_error)
}
