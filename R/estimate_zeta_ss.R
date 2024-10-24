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
}


#' Estimate \eqn{\zeta} using Smoothing Splines
#'
#' @param data A \code{list} or \code{data.table} object. Must contain variables/columns \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
#' @param df A \code{double} between \code{2} and \code{n}, where \code{n} is the number of unique elements of \code{SampleID}. If \code{df} is set to 2 or smaller, ordinary least squares regression is used to estimate \eqn{\zeta}, because this estimate is more reliable in these cases.
#' @param use_weights A \code{logical} value. If \code{TRUE}, weights are used to calculate \eqn{\zeta}. Weights are automatically determined.
#' @param simple_output A \code{logical} value. Set to \code{FALSE}, if all components of the zeta calculation is required.
#' @param na_rm A \code{logical} value. Should NA values be removed before estimating zeta hat?
#'
#' @return A \code{list} with one variable \code{zeta}, which is a double signifying the estimate of \eqn{\zeta}.
#' @export
#'
#' @examples print(1)

estimate_zeta_ss <- function(data, df = NULL, use_weights = FALSE, simple_output = TRUE, na_rm = TRUE){

  # Validate 'data' input
  if (!is.data.table(data)) {
    if (!is.list(data) && !is.data.frame(data)) {
      stop("Oops, it seems like 'data' is having an identity crisis. It dreams of being a data.table, list, or data.frame, but alas, it is currently trapped in the body of a sad and confused ", class(data), ". Poor data, always searching for its true purpose, like a GPS with no sense of direction!")
    }
    data <- as.data.table(data)
  }

  # Validate entries
  valdiate_estimate_zeta_ss(data = data, df = df, error = TRUE, funny = TRUE)

  # If df is 2, we use estimate_zeta() method
  if(!is.null(df) && df <= 2){
    return(estimate_zeta(data = data))
  }

  # Glocal precision estimates
  impr <- global_precision_estimates(data = data, silence = 1L)

  # Remove NA-values if one wishes to do so.
  if(na_rm){
    valid_ids <- (!is.na(data$MP_B)) & (!is.na(data$MP_A))
    data <- data[valid_ids, ]
  }

  # Calculate weight data
  weight_data <- NULL
  if(use_weights){
    R_i <- count_samplewise_replicates(data, summary = "none")$R_i
    vor_data <- fun_of_replicates2(data, "var") |> setDT()
    mor_data <- fun_of_replicates2(data, "mean") |> setDT()
    weight_data_mor <- weight_function(mor_data, vor_data, impr, output_type = "vector")
    weight_data <- sapply(1:length(mor_data$SampleID), function(x) rep(weight_data_mor[x], R_i[x]), simplify = FALSE) |> unlist()
  }
  else{
    weight_data <- rep(1, length(data$MP_A))
  }

  # Set order and assign variables based on lambda value
  axis_var <- if (impr$lambda < 0.5) {
    c("MP_A", "MP_B", "Var_A", "Var_B")
  } else {
    c("MP_B", "MP_A", "Var_B", "Var_A")
  }

  setorderv(data, axis_var[1], order = 1)
  x <- data[[axis_var[1]]]
  y <- data[[axis_var[2]]]
  var_h <- impr[[axis_var[3]]]
  var_v <- impr[[axis_var[4]]]
  x_unit <- x[order(x)]
  x_unit <- (x_unit - x_unit[1]) / diff(range(x_unit))
  all_knots <- c(0, calculate_interior_knots(x_unit), 1)

  # Smoothing spline fitting
  if(is.null(df)){
    ss_fit <- smooth.spline(x = x, y = y, w = weight_data, all.knots = all_knots, cv = TRUE, control.spar = list(low = 0.38, high = 1.5)) |> suppressWarnings()
  }
  else{
    ss_fit <- smooth.spline(x = x, y = y, w = weight_data, df = df, all.knots = all_knots, cv = FALSE)
  }

  var_eps <- ss_fit$pen.crit / (ss_fit$n - ss_fit$df)
  slopes <- predict(ss_fit, deriv = 1)$y
  var_h_trans <- mean(slopes^2, na.rm = TRUE) * var_h

  # Estimate zeta
  var_factor <- var_v + var_h_trans
  zeta_value <- var_eps * (ss_fit$n + ss_fit$df) / ss_fit$n / var_factor

  if(simple_output){
    return(list(zeta = zeta_value))
  }
  else{
    return(list(zeta = zeta_value, irr_var = var_factor, tot_var = var_eps, var_h = var_h, var_v = var_v, slopes = mean(slopes), sq_slopes = mean(slopes^2), shift = impr$lambda < 0.5))
  }
}

#' Estimate Comparison-wise \eqn{\zeta} using Smoothing Splines
#'
#' @param data A \code{data.table} object. Must contain variables/columns \code{comparison}, \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
#' @param df A \code{double} between \code{2} and \code{n}, where \code{n} is the number of unique elements of \code{SampleID}.
#'           Otherwise, it must be a \code{numeric} vector of length equal to the number of unique elements in \code{comparison}.
#'           If \code{df} is set to 2 or smaller, ordinary least squares regression is used to estimate \eqn{\zeta}, because this estimate is more reliable in these cases.
#' @param use_weights A \code{logical} value. If \code{TRUE}, weights are used to calculate \eqn{\zeta}. Weights are automatically determined.
#' @param na_rm A \code{logical} value. Should NA values be removed before estimating zeta hat?
#'
#' @return A \code{list} with one variable \code{zeta}, which is a double signifying the estimate of \eqn{\zeta}.
#' @export
#'
#' @examples print(1)

estimate_zetas_ss <- function(data, df = NULL, use_weights = TRUE, na_rm = TRUE){

  data_grouped_by_comparison <- split(data, by = "comparison", keep.by = TRUE)

  output_comparison_names <- names(data_grouped_by_comparison)
  output_comparison_ints <- 1L:length(data_grouped_by_comparison)


  if(!is.null(df)){
    if(length(unique(df)) == 1){
      df <- rep(df[1], length(data_grouped_by_comparison))
    }
    else if(!(length(df) == length(data_grouped_by_comparison))){
      stop("df must be of same length as the number of comparisons in data")
    }
  }

  output <- sapply(X = 1:length(data_grouped_by_comparison), FUN = function(x){
    estimate_zeta_ss(data = data_grouped_by_comparison[[x]],
                     df = df[x],
                     simple_output = TRUE,
                     na_rm = na_rm)
  }, simplify = FALSE)

  output <- rbindlist(output, idcol = "comparison")
  output$comparison <- output_comparison_names[match(output$comparison, output_comparison_ints)]

  return(output)

}

