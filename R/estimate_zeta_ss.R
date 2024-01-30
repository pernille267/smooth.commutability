#' Internal Function for validating entries
#'
#' This function is used internally validate entries in estimate_zeta_ss.
#' @keywords internal
valdiate_estimate_zeta_ss <- function(data, df, complex, advice, error = TRUE, funny = TRUE){

  # Check for required columns
  required_cols <- c("MP_A", "MP_B")
  if(!all(required_cols %in% names(data))){
    if(error){
      if(funny){
        stop("'data' must have variables named 'MP_B' and 'MP_A'.")
      }
      else{
        stop("'data' must have variables named 'MP_B' and 'MP_A'.")
      }
    }
    else{
      return(list(zeta = NA_real_))
    }
  }

  # Validation for 'df'
  if(!is.null(df) && (!is.numeric(df) || length(df) != 1)){
    if(error){
      if(funny){
        stop("Ah, the classic 'df' is at it again, trying to be something it's not! 'df' should be a non-missing numeric value or NULL, but it's off gallivanting as a ", class(df), ". It's like a wannabe actor auditioning for a role it was never meant to play!")
      }
      else{
        stop("'df' is expected to be a non-missing numeric value or NULL, but is a ", class(df), ".")
      }
    }
    else{
      return(list(zeta = NA_real_))
    }
  }

  # Validation for 'complex'
  if (!is.logical(complex) || length(complex) != 1){
    if(error){
      if(funny){
        stop("'complex' must be a non-missing logical value (TRUE/FALSE), but it is tragically wearing the guise of a ", class(complex), ". The universe can be so cruel and complex sometimes, can't it?")
      }
      else{
        stop("'complex' must be a non-missing logical value (TRUE/FALSE), but is a ", class(complex), ".")
      }
    }
    else{
      return(list(zeta = NA_real_))
    }
  }

  # Validation for 'advice'
  if (!is.logical(advice) || length(advice) != 1){
    if(error){
      if(funny){
        stop("'advice' must be a non-missing logical value (TRUE/FALSE). Remember, computers are like Old Testament gods; lots of rules and no mercy.")
      }
      else{
        stop("'advice' must be a non-missing logical value (TRUE/FALSE), but is a ", class(advice), ".")
      }
    }
    else{
      return(list(zeta = NA_real_))
    }
  }
}


#' Estimate \eqn{\zeta} using Smoothing Splines
#'
#' @param data A \code{list} or \code{data.table} object. Must contain variables/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
#' @param df A \code{double} between \code{2} and \code{n}, where \code{n} is the number of unique elements of \code{SampleID}. If \code{df} is set to 2 or smaller, ordinary least squares regression is used to estimate \eqn{\zeta}, because this estimate is more reliable in these cases.
#' @param complex A \code{logical} value. Should the Delta method be used to approximate the variance of the f? If yes, set to \code{TRUE}.
#' @param advice A \code{logical} value. Should a friendly piece of advice be delivered to you? Set to \code{TRUE} if you want this.
#'
#' @return A \code{list} with one variable \code{zeta}, which is a double signifying the estimate of \eqn{\zeta}.
#' @export
#'
#' @examples print(1)

estimate_zeta_ss <- function(data, df = NULL, complex = TRUE, advice = FALSE){

  # Validate 'data' input
  if (!is.data.table(data)) {
    if (!is.list(data) && !is.data.frame(data)) {
      stop("Oops, it seems like 'data' is having an identity crisis. It dreams of being a data.table, list, or data.frame, but alas, it is currently trapped in the body of a sad and confused ", class(data), ". Poor data, always searching for its true purpose, like a GPS with no sense of direction!")
    }
    data <- as.data.table(data)
  }

  # Validate entries
  valdiate_estimate_zeta_ss(data = data, df = df, complex = complex, advice = advice, error = TRUE, funny = TRUE)

  if(!is.null(df) && df <= 2){
    return(estimate_zeta(data = data))
  }

  # Glocal precision estimates
  impr <- global_precision_estimates(data = data, silence = 1L)

  # Set order and assign variables based on lambda value
  axis_var <- if (impr$lambda < 0.9) {
    if (advice) cat("Advice: Ratio of measurement error variances is < 0.9. Consider switching axes.\n")
    c("MP_A", "MP_B", "Var_A", "Var_B")
  } else {
    c("MP_B", "MP_A", "Var_B", "Var_A")
  }

  setorderv(data, axis_var[1], order = 1)
  x <- data[[axis_var[1]]]
  y <- data[[axis_var[2]]]
  var_h <- impr[[axis_var[3]]]
  var_v <- impr[[axis_var[4]]]

  # Smoothing spline fitting
  if(is.null(df)){
    cv_crit <- function(df, x, y){
      smooth.spline(x = x, y = y, df = df, all.knots = TRUE)$cv.crit
    }
    df <- optimise(f = cv_crit, x = x, y = y, interval = c(2, length(x)))$minimum
    df <- min(df, 10)
    ss_fit <- smooth.spline(x = x, y = y, df = df, all.knots = TRUE)
  }
  else{
    ss_fit <- smooth.spline(x = x, y = y, df = df, all.knots = TRUE)
  }

  var_eps <- ss_fit$pen.crit / (ss_fit$n - ss_fit$df)
  var_h_trans <- mean(predict(ss_fit, deriv = 1)$y^2) * var_h

  # Estimate zeta
  var_factor <- if (complex) var_v + var_h_trans else var_v + var_h
  zeta_value <- var_eps * (ss_fit$n + ss_fit$df) / ss_fit$n / var_factor

  return(list(zeta = zeta_value))

}
