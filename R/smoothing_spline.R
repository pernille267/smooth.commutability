#' Internal Function for Data Conversion
#'
#' This function is used internally for converting different data types to data.table.
#' @keywords internal
convert_to_data_table <- function(data){
  if(is.list(data)){
    setDT(data)
  }
  else if(is.data.frame(data)){
    as.data.table(data)
  }
  else{
    stop("'data' is expected to be a data.table object but is a ", class(data))
  }
}

#' Calculate interior knots
#'
#' This function is used internally for calculating interior knots.
#' @keywords internal
calculate_interior_knots <- function(x) {
  unique_x <- unique(x)
  n_interior_knots <- length(unique_x)
  interior_knots <- quantile(x = unique_x, probs = seq(from = 0, to = 1, length.out = n_interior_knots + 2))[-c(1, n_interior_knots + 2)]
  names(interior_knots) <- NULL
  return(interior_knots)
}

#' Fit a Smoothing Spline Model with O'Sullivan Penalty
#'
#' This function fits a smoothing spline model to the provided data using an O'Sullivan penalty matrix. It is designed to work with various data types and offers flexibility in model tuning through degrees of freedom (`df`) and penalty parameter (`lambda`).
#'
#' @param data A \code{data.table}, \code{data.frame} or \code{list} containing the measurements.
#'             It should include two specific variables named 'MP_B' and 'MP_A'.
#' @param df   Optional; a \code{numeric} value specifying the degrees of freedom for the spline model.
#'             Should be between 2 and the number of clinical samples in \code{data}.
#'             The default value (\code{NULL}) triggers automatic selection using n-fold cross-validation.
#' @param lambda Optional; a non-negative \code{numeric} value specifying the penalty strength.
#'               The default (\code{NULL}) results in the function choosing an optimal value based on
#'               n-fold cross-validation to minimize the mean squared prediction error.
#' @param df_max A \code{numeric} value between 2 and the number of clinical samples in \code{data}.
#'               If n-fold cross-validation is used to obtain the optimal effective number of degrees of freedom, what is its upper limit.
#'
#' @return A \code{list} containing the components of the smoothing spline fit, including:
#'         - \code{fitted}: Fitted values of the spline model.
#'         - \code{residuals}: Residuals from the fitted model.
#'         - \code{df}: Effective degrees of freedom used.
#'         - \code{lambda}: Penalty parameter used.
#'         - \code{cv_crit}: Criterion value from cross-validation.
#'         - Additional model components.
#'         These components are often utilized in subsequent analyses and visualizations.
#'
#' @export
#'
#' @examples
#' # Example using a dataset with columns MP_B and MP_A
#' library(data.table)
#' data_example <- data.table(MP_B = 1:10, MP_A = rnorm(10))
#' result <- smoothing_spline(data_example, df = 5)
#' print(result$fitted)

smoothing_spline <- function(data, df = NULL, lambda = NULL, df_max = 5){

  MP_B <- NULL
  cv_crit <- NA
  data <- convert_to_data_table(data)

  # Check for required variables in data
  if(!any("MP_B" == names(data)) | !any("MP_A" == names(data))){
    stop("'data' is expected to have variables named 'MP_B' and 'MP_A'.")
  }

  # Sort according to predictor values
  setorder(x = data, MP_B)
  orig_x <- data$MP_B
  y <- data$MP_A

  # Transform predictor values to [0, 1] to avoid overflow
  x <- (orig_x - orig_x[1]) / diff(range(orig_x))

  # Calculation of interior knots and B-spline basis matrix
  interior_knots <- calculate_interior_knots(x)
  x_min_max <- range(x)
  B <- bs(x = x, knots = interior_knots, degree = 3L, intercept = TRUE, Boundary.knots = x_min_max)

  # Penalty matrix and cross-products
  Omega <- penalty_matrix(x_min = x_min_max[1], x_max = x_min_max[2], interior_knots = interior_knots)
  BTB <- crossprod(B, B)
  BTy <- crossprod(B, y)

  # Degrees of freedom / smoothing parameter calculations

  # both "df" and "lambda" is missing
  if(is.null(df) && is.null(lambda)){
    # Obtain the "spar" that minimizes the n-fold cross-validation curve
    nfold_cv_opt <- optimise(f = nfold_cv, lower = -1.5, upper = 1.5, y = y, B = B, BTB = BTB, Omega = Omega, cross_validation = TRUE)
    # Map "spar" to the smoothing parameter "lambda"
    r <- sum(diag(BTB)) / sum(diag(Omega))
    lambda <- r * 256**(3 * nfold_cv_opt$minimum - 1)
    # If df0 exceeds the maximum allowable limit (i.e., df_max), obtain a new "lambda" based "df_max".
    df0 <- sum(diag(B %*% solve(BTB + lambda * Omega) %*% t(B)))
    if(df0 > df_max){
      lambda <- to_lambda(df = df_max, B = B, BTB = BTB, Omega = Omega)
      cv_crit <- nfold_cv(lambda = lambda, y = y, B = B, BTB = BTB, Omega = Omega)
    }
    else{
      cv_crit <- nfold_cv_opt$objective
    }
  }

  # "lambda" is missing, but "df" exists
  else if(!is.null(df) && is.null(lambda)){
    lambda <- to_lambda(df = df, B = B, BTB = BTB, Omega = Omega)
    cv_crit <- nfold_cv(lambda = lambda, y = y, B = B, BTB = BTB, Omega = Omega)
  }

  # "df" is missing, but "lambda" exists
  else if(is.null(df) && !is.null(lambda)){
    cv_crit <- nfold_cv(lambda = lambda, y = y, B = B, BTB = BTB, Omega = Omega)
  }

  # Spline fitting
  Q <- solve(BTB + lambda * Omega)
  beta <- Q %*% BTy
  S <- B %*% Q %*% t(B)
  df <- sum(diag(S))
  fitted <- S %*% y
  residuals <- y - fitted
  var_eps <- sum(residuals^2) / (length(x) - df)
  cov_beta <- var_eps * Q %*% t(B) %*% B %*% Q

  # Returning the fitted object as a list
  return(list(fitted = fitted,
              residuals = residuals,
              df = df,
              lambda = lambda,
              cv_crit = cv_crit,
              var_eps = var_eps,
              coefficients = beta,
              cov_beta = cov_beta,
              smoothing_matrix = S,
              penalty_matrix = Omega,
              B = B,
              interior_knots = interior_knots,
              knot_boundary = x_min_max))
}
