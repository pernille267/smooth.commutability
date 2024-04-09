#' Prediction Interval Estimation for External Quality Assessment Data via Smoothing Splines
#'
#' @param data A \code{data.table}, \code{list} or \code{data.frame} object. Should include the ID column \code{SampleID} and the measurement columns \code{MP_A} and \code{MP_B}. Clinical sample measurements should be in here.
#' @param new_data A \code{data.table}, \code{list} or \code{data.frame} object. Can include the ID column \code{SampleID} and the measurement column \code{MP_A}, but the measurement column \code{MP_B} is mandatory. External quality assessment (EQA) material measurements should be in here.
#' @param df A optional \code{numeric} value greater than or equal to 2, but less than or equal to the number of clinical samples in \code{data}. If both \code{df} and \code{lambda} is set to \code{NULL}, n-fold cross-validation is used to obtain the optimal \code{df}.
#' @param lambda A optional \code{numeric} value greater than 0. If both \code{df} and \code{lambda} is set to \code{NULL}, n-fold cross-validation is used to obtain the optimal \code{lambda}.
#' @param df_max A \code{numeric} value between 2 and the number of clinical samples in \code{data}. If n-fold cross-validation is used to obtain the optimal effective number of degrees of freedom, what is its upper limit.
#' @param R_ratio A \code{numeric} value indicating the ratio of replicates between that number \code{data} and that number \code{new_data} are based on. Only relevant if the number of replicates of \code{data} and \code{new_data} differs.
#' @param level A \code{numeric} value representing the confidence level for the approximated prediction intervals. It should be between \code{0} and \code{1}. The default setting is \code{0.99}. Please adjust for simultaneous testing if pointwise prediction intervals are used for classifying more than one EQA material / reference material in the same IVD-MD comparison.
#' @param rounding An \code{integer} specifying the desired decimal places for the predictions and prediction intervals. The default setting is three, offering sufficient precision. The maximum limit is twelve due to the utilization of double numbers.
#'
#' @details For optimal results, we recommended to include all \code{SampleID}, \code{MP_A}, and \code{MP_B} in \code{new_data} whenever possible. When constructing prediction bands, only \code{MP_B} should be available.
#'
#' @return A \code{data.table} object comprising approximated prediction interval data based on user inputs.
#' @export
#'
#' @examples print(1)

predict_smoothing_spline <- function(data, new_data, df = NULL, lambda = NULL, df_max = 7.5, R_ratio = 1, level = 0.99, rounding = 3L){

  MP_B <- NULL

  # Validate 'data' input
  #-----------------------------------------------------------------------------
  if (!is.data.table(new_data)) {
    if (!is.list(new_data) && !is.data.frame(new_data)) {
      stop("Oops, it seems like 'new_data' is having an identity crisis. It dreams of being a data.table, list, or data.frame, but alas, it is currently trapped in the body of a sad and confused ", class(new_data), ". Poor 'new_data', always searching for its true purpose, like a GPS with no sense of direction!")
    }
    new_data <- as.data.table(new_data)
  }

  if(!any(names(new_data) == "MP_B")){
    stop("MP_B is missing in 'new_data'. Make sure it is part of 'new_data'.")
  }

  # Sort according to new predictor values
  new_data_clone <- copy(new_data)
  old_order <- order(new_data_clone$MP_B)
  setorder(x = new_data_clone, MP_B)
  orig_nx <- new_data_clone$MP_B

  # Get relevant components from smoothing spline fit
  smoothing_spline_fit <- smoothing_spline(data = data, df = df, lambda = lambda, df_max = df_max)
  var_eps <- smoothing_spline_fit$var_eps
  cov_beta <- smoothing_spline_fit$cov_beta
  interior_knots <- smoothing_spline_fit$interior_knots
  x_min_max <- smoothing_spline_fit$knot_boundary
  df <- smoothing_spline_fit$df
  beta <- smoothing_spline_fit$coefficients


  nx <- (orig_nx - min(data$MP_B)) / diff(range(data$MP_B))
  outside_id <- integer()
  if(all(nx < 0 | nx > 1)){
    stop("Commutability evaluation outside the concentration interval of clinical samples is currently not possible.")
  }
  else if(any(nx < 0 | nx > 1)){
    warning("Commutability evaluation outside the concentration interval of clinical samples is currently not possible.")
    outside_id <- which(nx < 0 | nx > 1)
    nx <- nx[-outside_id]
    new_data_clone <- new_data_clone[-outside_id, ]
  }

  # Get variance of new predictor values
  B_new <- bs(x = nx, knots = interior_knots, degree = 3L, intercept = TRUE, Boundary.knots = x_min_max)
  var_new_values <- diag(B_new %*% cov_beta %*% t(B_new))

  # Get variance of prediction error
  var_pred_error <- (var_eps + var_new_values) * R_ratio

  # Ouptut
  n <- length(data$MP_B[!is.na(data$MP_B)])
  t_quant <- qt(p = (1 - level) / 2, df = n - df, lower.tail = FALSE)
  ny <- round((B_new %*% beta)[,], rounding)
  lwr <- round(pmax(0, ny - t_quant * sqrt(var_pred_error)), rounding)
  upr <- round(ny + t_quant * sqrt(var_pred_error), rounding)
  nx <- round(nx, rounding)

  if(all(c("SampleID", "MP_A") %in% names(new_data_clone))){
    output <- list("SampleID" = new_data_clone$SampleID, "MP_B" = round(new_data_clone$MP_B, rounding), "MP_A" = round(new_data_clone$MP_A, rounding), "prediction" = ny, "lwr" = lwr, "upr" = upr)
    output$inside <- as.integer(output$MP_A >= output$lwr & output$MP_A <= output$upr)
  }
  else{
    if("SampleID" %in% names(new_data_clone)){
      output <- list("SampleID" = new_data_clone$SampleID, "MP_B" = nx, "prediction" = ny, "lwr" = lwr, "upr" = upr)
    }
    else if("MP_A" %in% names(new_data_clone)){
      output <- list("MP_B" = round(new_data_clone$MP_B, rounding), "MP_A" = round(new_data_clone$MP_A, rounding), "prediction" = ny, "lwr" = lwr, "upr" = upr)
      output$inside <- as.integer(output$MP_A >= output$lwr & output$MP_A <= output$upr)
    }
    else{
      output <- list("MP_B" = round(new_data_clone$MP_B, rounding), "prediction" = ny, "lwr" = lwr, "upr" = upr)
    }
  }

  setDT(output)
  output <- output[old_order,]
  return(output)

}
