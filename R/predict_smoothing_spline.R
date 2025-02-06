#' Check input of predict_smoothing_spline
#'
#' This function is used internally for error handling regarding bad input.
#' @keywords internal
validate_prediction_data_components <- function(data, new_data){

  # Defaults
  convert_data <- FALSE
  required_columns_data <- c("MP_B", "MP_A", "SampleID")
  convert_new_data <- FALSE

  # Validate 'data' argument
  #-----------------------------------------------------------------------------

  # Check if class if class is correct
  if (!is.data.table(data)) {
    if (!is.list(data) && !is.data.frame(data)) {
      stop("The data argument is not a data.table, data.frame or list, but a: ", class(data))
    }
    convert_data <- TRUE
  }
  missing_columns_data <- required_columns_data[which(!(required_columns_data %in% names(data)))]

  # Check if have the correct names
  if(length(missing_columns_data) >= 1){
    stop("The data argument does not have the necessary columns. ", paste(missing_columns_data, collapse = ", "), " are missing.")
  }

  # Validate 'new_data' argument
  #-----------------------------------------------------------------------------

  # Check if class if class is correct
  if (!is.data.table(new_data)) {
    if (!is.list(new_data) && !is.data.frame(new_data)) {
      stop("The new_data argument is not a data.table, data.frame or list, but a: ", class(new_data))
    }
    convert_new_data <- TRUE
  }

  # Check if have the correct names
  if(!any("MP_B" == names(new_data))){
    stop("The new_data argument does not have the necessary columns. MP_B is missing.")
  }

  return(c(convert_data, convert_new_data))

}

#' Prediction Interval Estimation for External Quality Assessment Data via Smoothing Splines
#'
#' @param data A \code{data.table}, \code{list} or \code{data.frame} object.
#'             Should include the ID column \code{SampleID} and the measurement
#'             columns \code{MP_A} and \code{MP_B}.
#'             Clinical sample measurements should be in here.
#' @param new_data A \code{data.table}, \code{list} or \code{data.frame} object.
#'                 Can include the ID column \code{SampleID} and the measurement
#'                 column \code{MP_A}, but the measurement column \code{MP_B} is mandatory.
#'                 External quality assessment (EQA) material measurements should be in here.
#' @param weighted A \code{logical} value. If \code{TRUE}, weights are
#'                 iteratively estimated from data.
#'                 If set to \code{NULL}, no weights are used. See details for more information.
#' @param df A optional \code{numeric} value greater than or equal to 2,
#'           but less than or equal to the number of clinical samples in \code{data}.
#'           If both \code{df} and \code{lambda} is set to \code{NULL},
#'           n-fold cross-validation is used to obtain the optimal \code{df}.
#' @param lambda A optional \code{numeric} value greater than 0.
#'               If both \code{df} and \code{lambda} is set to \code{NULL},
#'               n-fold cross-validation is used to obtain the optimal \code{lambda}.
#' @param df_max A \code{numeric} value between 2 and the number of clinical samples in \code{data}.
#'               If n-fold cross-validation is used to obtain an optimal
#'               effective number of degrees of freedom. This value is
#'               capped at this upper limit.
#' @param R_ratio A \code{numeric} value indicating the ratio of replicates between
#'                that number \code{data} and that number \code{new_data} are based on.
#'                Only relevant if the number of replicates of \code{data} and \code{new_data} differs.
#' @param level A \code{numeric} value representing the confidence level for the
#'              approximated prediction intervals. It should be between \code{0} and \code{1}.
#'              The default setting is \code{0.99}.
#' @param simultaneous A \code{logical} value. If set to \code{TRUE},
#'                     prediction intervals are Bonferroni-corrected.
#'                     The default is set to \code{FALSE}.
#' @param negative_ok A \code{logical} value. If set to \code{TRUE}, negative values are allowed.
#'                    See details.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#' @param include_prediction_variance A \code{logical} value. If \code{TRUE}, prediction variances are returned part of the output.
#' @param rounding An \code{integer} specifying the desired decimal places for the
#'                 predictions and prediction intervals. The default setting is \code{3L},
#'                 offering sufficient precision. The maximum limit is twelve
#'                 due to double precision.
#' @param additional_parameters A \code{list} containing additional alternatives for the smoothing spline fit.
#'
#' @description
#' Predict and estimate prediction intervals based on new observations using
#' smoothing splines.
#'
#' @details
#' Structure
#'
#' The input arguments, \code{data}, \code{new_data} and \code{weight_data} have
#' rather strict requirements for them to be accepted as valid inputs. Here are
#' the requirements for each of these three arguments.
#'
#' The argument \code{data} must be a \code{data.table}, \code{list} or \code{data.frame}
#' with the following columns:
#' \itemize{
#'    \item \code{SampleID}: The sample identifiers for each sample in the dataset
#'    \item \code{MP_A}: The measurements of the first IVD-MD in the IVD-MD comparison. Response variable.
#'    \item \code{MP_B}: The measurements of the second IVD-MD in the IVD-MD comparison. Predictor variable.
#' }
#'
#' The argument \code{new_data} must be a \code{data.table}, \code{list} or \code{data.frame}
#' with the following columns:
#'
#' \itemize{
#'    \item \code{SampleID}*: The sample identifiers for each sample in the dataset
#'    \item \code{MP_A}*: The measurements of the first IVD-MD in the IVD-MD comparison. Actual response variable.
#'    \item \code{MP_B}: The measurements of the second IVD-MD in the IVD-MD comparison. New predictor variable.
#' }
#' The columns \code{SampleID} and \code{MP_A} are optional, but required for
#' calculation of \code{inside}.
#'
#' The argument \code{weighted} can be either \code{TRUE} or \code{FALSE}.
#' If \code{FALSE}, no weights are used in prediction.
#' If \code{TRUE}, weights are iteratively estimated based on data.
#'
#' Smoothing parameters
#'
#' The arguments \code{lambda}, \code{df} and \code{df_max} are meant to control
#' the smoothness of the fitted smoothing spline. The two former, \code{lambda}
#' and \code{df} are straight-forward. Here \code{lambda} is the computational
#' smooting parameter and \code{df} is the effective number of degrees of freedom.
#' There are four possibilities of combinations of \code{lambda} and \code{df}
#' \itemize{
#'    \item Both \code{lambda} and \code{df} are \code{NULL}. Then leave-one-out Cross-validation
#'    is used to obtain an optimal value of \code{lambda}, and derive the corresponding
#'    value of \code{df}
#'    \item \code{lambda} is \code{NULL}, but \code{df} is given. Then, numerical optimization
#'    is used to obtain the \code{lambda} that corresponds with the given \code{df}.
#'    \item \code{lambda} is given, but \code{df} is \code{NULL}. Then no optimization
#'    is necessary.
#'    \item Both \code{lambda} and \code{df} are given. \code{df} is then ignored
#'    and \code{lambda} is used.
#' }
#' If Both \code{lambda} and \code{df} are \code{NULL}. Then, \code{df_max} is
#' used to control how large the optimal value of \code{df} can be based on
#' leave-one-out Cross-validation. If the optimal value is larger than \code{df_max},
#' this value is used in favor of the optimal. In such cases, numerical optimization
#' is used to obtain the \code{lambda} that corresponds with the given \code{df_max}.
#'
#' Prediction interval estimation
#'
#' For the prediction interval estimation, the width of the estimated prediction
#' intervals can be modified using \code{level}, \code{simultaneous} and \code{R_ratio}
#' arguments.
#'
#' The former, \code{level}, controls the required confidence level
#' of the estimated prediction intervals. This must of course be a value between \code{0}
#' and \code{1}, but typical values are \code{0.90}, \code{0.95} and \code{0.99}.
#'
#' The \code{simultaneous} argument controls whether confidence levels should be
#' Bonferroni-corrected for when estimating multiple prediction intervals at once. The
#' effective point-wise confidence level for a particular value of \code{level}
#' is given by \eqn{\mathrm{level}_{\mathrm{eff}} = 1 - (1 - \mathrm{level}) / m},
#' where \code{m} is the number of joint prediction intervals.
#'
#' The \code{R_ratio} argument controls the ratio of the number of replicates used
#' in \code{data} compared to what is used in \code{new_data}. If more replicates
#' are used in \code{new_data} compared to \code{data}, the widths of the estimated
#' prediction interval should be narrower.
#'
#' Negative results
#'
#' For standard measurement results, we do not expect them to be negative.
#' Nevertheless, sometimes, the predictions and estimated prediction intervals
#' may be negative. If you wish to cutoff negative results at \code{zero}, you
#' can do this by setting \code{negative_ok} to \code{FALSE}.
#'
#' However, in some cases, you might expect negative results. For example,
#' when data is log-transformed and some of the raw measurements are between
#' \code{0} and \code{1}. In these cases, it might be favorable to set \code{negative_ok} to \code{TRUE}.
#'
#'
#' @return A \code{data.table} object comprising approximated prediction interval data based on user inputs.
#' @export
#'
#' @examples print(1)

predict_smoothing_spline <- function(data,
                                     new_data,
                                     weighted = FALSE,
                                     df = NULL,
                                     lambda = NULL,
                                     df_max = 7.5,
                                     R_ratio = 1,
                                     level = 0.99,
                                     simultaneous = FALSE,
                                     negative_ok = TRUE,
                                     attempt_fast = FALSE,
                                     include_prediction_variance = FALSE,
                                     rounding = 3L,
                                     additional_parameters = list("method" = "loocv",
                                                                  "smudge" = 1.4,
                                                                  "c_star_min" = 0.38,
                                                                  "c_star_max" = 1.50,
                                                                  "tol" = 0.5,
                                                                  "window" = 10,
                                                                  "iter" = 1)){

  # Defaults
  MP_B <- NULL
  weights <- NULL

  # Check validity of data, new_data and weight_data
  validate_and_fix_input <- validate_prediction_data_components(data, new_data)

  # Conversions if necessary and weight calculations
  if(validate_and_fix_input[1]){
    data <- as.data.table(data)
  }
  if(validate_and_fix_input[2]){
    new_data <- as.data.table(new_data)
  }

  # Handling training weighting..
  if(weighted){
    weights <- "estimate"
  }
  else{
    weights <- rep(1, length(data$MP_B))
  }

  # Sort according to new predictor values
  new_data_clone <- copy(new_data)
  old_order <- order(new_data_clone$MP_B, na.last = FALSE)
  setorder(x = new_data_clone, MP_B, na.last = FALSE)
  nx <- new_data_clone$MP_B

  old_order <- match(order(nx, na.last = FALSE), old_order)

  # Get relevant components from smoothing spline fit
  smoothing_spline_fit <- smoothing_spline(data = data,
                                           weights = weights,
                                           df = df,
                                           lambda = lambda,
                                           df_max = df_max,
                                           attempt_fast = attempt_fast,
                                           additional_parameters = additional_parameters)
  var_eps <- smoothing_spline_fit$var_eps
  cov_beta <- smoothing_spline_fit$cov_beta
  interior_knots <- smoothing_spline_fit$interior_knots
  x_min_max <- smoothing_spline_fit$knot_boundary
  df <- smoothing_spline_fit$df
  beta <- smoothing_spline_fit$coefficients
  n <- length(data$MP_B[!is.na(data$MP_B)])
  m <- if(isTRUE(simultaneous)){length(new_data$MP_B[!is.na(new_data$MP_B)])}else{1}

  # Compute interpolation and extrapolation predictions and their variances
  nu <- (nx - min(data$MP_B, na.rm = TRUE)) / diff(range(data$MP_B, na.rm = TRUE))
  nu_na <- any(is.na(nu))
  na_id <- NULL
  if(isTRUE(nu_na)){
    na_id <- which(is.na(nu))
    nu <- nu[-na_id]
  }

  extrapolate_lower <- nu < 0
  extrapolate_upper <- nu > 1
  extrapolate <- extrapolate_lower | extrapolate_upper
  interpolate <- !(extrapolate)
  ny <- nu
  nw <- nu
  var_new_values <- nu
  extrapolation <- nu

  # Handling new weighting..
  if(weighted){
    nw_model <- smooth.spline(x = smoothing_spline_fit$u,
                              y = log(smoothing_spline_fit$weights),
                              df = n/2)
    nw <- exp(predict(nw_model, x = nu)$y)
  }
  else{
    nw <- rep(1, length(nu))
  }

  if(any(interpolate)){
    B_new <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)),
                        x = nu[interpolate],
                        derivs = rep(0, length(nu[interpolate])),
                        outer.ok = TRUE)$design
    ny[interpolate] <- round((B_new %*% beta)[,], rounding)
    var_new_values[interpolate] <- calculate_pred_var(B_new, cov_beta)[,]
    extrapolation[interpolate] <- 0
  }

  if(any(extrapolate)){
    extrapolation[extrapolate] <- 1
    B_boundary <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)),
                             x = c(0, 1),
                             derivs = c(0, 0),
                             outer.ok = TRUE)$design
    B_boundary_deriv <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)),
                                   x = c(0, 1),
                                   derivs = c(1, 1),
                                   outer.ok = TRUE)$design
    boundary_knot_derivs <- B_boundary_deriv %*% beta |> as.vector()
    boundary_knot_fitted <- range(smoothing_spline_fit$fitted)
    boundary_knot_derivs_vars <- calculate_pred_var(B_boundary_deriv, cov_beta)[,]
    boundary_knot_fitted_vars <- calculate_pred_var(B_boundary, cov_beta)[,]
    if(any(extrapolate_lower)){
      ny[extrapolate_lower] <- round(boundary_knot_fitted[1L] + boundary_knot_derivs[1L] * (nu[extrapolate_lower] - 0), rounding)
      var_new_values[extrapolate_lower] <- boundary_knot_fitted_vars[1L] + (nu[extrapolate_lower] - 0)^2 * boundary_knot_derivs_vars[1L]
    }
    if(any(extrapolate_upper)){
      ny[extrapolate_upper] <- round(boundary_knot_fitted[2L] + boundary_knot_derivs[2L] * (nu[extrapolate_upper] - 1), rounding)
      var_new_values[extrapolate_upper] <- boundary_knot_fitted_vars[2L] + (nu[extrapolate_upper] - 1)^2 * boundary_knot_derivs_vars[2L]
    }
  }

  # Get variance of prediction error
  var_pred_error <- (var_eps / nw + var_new_values) * R_ratio

  # Ouptut
  negative_results <- any(data$MP_A < 0, na.rm = TRUE) | any(data$MP_B < 0, na.rm = TRUE)
  t_quant <- qt(p = (1 - level) / 2 / m, df = n - df, lower.tail = FALSE)
  if(isFALSE(negative_ok)){
    lwr <- round(pmax(c(-1e6, 0)[c(negative_results, !negative_results)], ny - t_quant * sqrt(var_pred_error)), rounding)
    upr <- round(pmax(c(-1e6, 0)[c(negative_results, !negative_results)], ny + t_quant * sqrt(var_pred_error)), rounding)
  }
  else{
    lwr <- round(ny - t_quant * sqrt(var_pred_error), rounding)
    upr <- round(ny + t_quant * sqrt(var_pred_error), rounding)
  }

  nu <- round(nu, rounding)

  if(all(c("SampleID", "MP_A") %in% names(new_data_clone))){

    if(isTRUE(nu_na)){
      ny <- c(rep(NA, length(na_id)), ny)
      lwr <- c(rep(NA, length(na_id)), lwr)
      upr <- c(rep(NA, length(na_id)), upr)
      extrapolation <- c(rep(NA, length(na_id)), extrapolation)
      var_pred_error <- c(rep(NA, length(na_id)), var_pred_error)
    }

    output <- list("SampleID" = new_data_clone$SampleID,
                   "MP_B" = round(new_data_clone$MP_B, rounding),
                   "MP_A" = round(new_data_clone$MP_A, rounding),
                   "prediction" = ny,
                   "lwr" = lwr,
                   "upr" = upr)

    output$inside <- as.integer(output$MP_A >= output$lwr & output$MP_A <= output$upr)
    output$extrapolate <- extrapolation
    output$var_pred_error <- if(include_prediction_variance){var_pred_error}else{NULL}
  }
  else{
    if("SampleID" %in% names(new_data_clone)){

      if(isTRUE(nu_na)){
        ny <- c(rep(NA, length(na_id)), ny)
        lwr <- c(rep(NA, length(na_id)), lwr)
        upr <- c(rep(NA, length(na_id)), upr)
        extrapolation <- c(rep(NA, length(na_id)), extrapolation)
        var_pred_error <- c(rep(NA, length(na_id)), var_pred_error)
      }

      output <- list("SampleID" = new_data_clone$SampleID,
                     "MP_B" = nx,
                     "prediction" = ny,
                     "lwr" = lwr,
                     "upr" = upr)
      output$extrapolate <- extrapolation
      output$var_pred_error <- if(include_prediction_variance){var_pred_error}else{NULL}
    }
    else if("MP_A" %in% names(new_data_clone)){

      if(isTRUE(nu_na)){
        ny <- c(rep(NA, length(na_id)), ny)
        lwr <- c(rep(NA, length(na_id)), lwr)
        upr <- c(rep(NA, length(na_id)), upr)
        extrapolation <- c(rep(NA, length(na_id)), extrapolation)
        var_pred_error <- c(rep(NA, length(na_id)), var_pred_error)
      }

      output <- list("MP_B" = round(new_data_clone$MP_B, rounding),
                     "MP_A" = round(new_data_clone$MP_A, rounding),
                     "prediction" = ny,
                     "lwr" = lwr,
                     "upr" = upr)

      output$inside <- as.integer(output$MP_A >= output$lwr & output$MP_A <= output$upr)
      output$extrapolate <- extrapolation
      output$var_pred_error <- if(include_prediction_variance){var_pred_error}else{NULL}
    }
    else{

      if(isTRUE(nu_na)){
        ny <- c(rep(NA, length(na_id)), ny)
        lwr <- c(rep(NA, length(na_id)), lwr)
        upr <- c(rep(NA, length(na_id)), upr)
        extrapolation <- c(rep(NA, length(na_id)), extrapolation)
        var_pred_error <- c(rep(NA, length(na_id)), var_pred_error)
      }

      output <- list("MP_B" = round(new_data_clone$MP_B, rounding),
                     "prediction" = ny,
                     "lwr" = lwr,
                     "upr" = upr)
      output$extrapolate <- extrapolation
      output$var_pred_error <- if(include_prediction_variance){var_pred_error}else{NULL}
    }
  }

  setDT(output)

  output <- output[old_order,]
  #class(output) <- "predict_smoothing_spline"
  return(output)
}

#' Comparison-wise Prediction Interval Estimation for External Quality Assessment Data via Smoothing Splines
#'
#' @param data A grouped \code{data.table}, \code{list} or \code{data.frame}
#'             object by \code{comparison}. Should include the ID column \code{SampleID}
#'             and the measurement columns \code{MP_A} and \code{MP_B}.
#'             Clinical sample measurements should be in here.
#' @param new_data A grouped \code{data.table}, \code{list} or \code{data.frame}
#'                 object by \code{comparison}. Can include the ID column \code{SampleID}
#'                 and the measurement column \code{MP_A}, but the measurement column \code{MP_B} is mandatory.
#'                 External quality assessment (EQA) material measurements should be in here.
#' @param weighted A \code{logical} value. If \code{TRUE}, weights are
#'                 iteratively estimated from data.
#'                 If set to \code{NULL}, no weights are used. See details for more information.
#' @param df A optional \code{numeric} vector of length equal to the number of
#'           unique IVD-MD comparisons in \code{data}. It can alternatively be
#'           a single \code{double} value or \code{NULL}, which then applies to all IVD-MD comparions.
#' @param lambda A optional \code{numeric} vector of length equal to the number of
#'               unique IVD-MD comparisons in \code{data}. It can alternatively be
#'               a single \code{double} value or \code{NULL}, which then applies to all IVD-MD comparions.
#' @param df_max A \code{numeric} value between 2 and the number of clinical samples in \code{data}.
#' @param R_ratio A \code{numeric} value indicating the ratio of replicates between
#'                that number \code{data} and that number \code{new_data} are based on.
#'                Only relevant if the number of replicates of \code{data} and \code{new_data} differs.
#' @param level A \code{numeric} value representing the confidence level for the
#'              approximated prediction intervals. It should be between \code{0} and \code{1}.
#'              The default setting is \code{0.99}.
#' @param simultaneous A \code{logical} value. If set to \code{TRUE},
#'                     prediction intervals are Bonferroni-corrected.
#'                     The default is set to \code{FALSE}.
#' @param negative_ok A \code{logical} value. If set to \code{TRUE}, negative values are allowed.
#'                    See details.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#' @param include_prediction_variance A \code{logical} value. If \code{TRUE}, prediction variances are returned part of the output.
#' @param rounding An \code{integer} specifying the desired decimal places for the
#'                 predictions and prediction intervals. The default setting is \code{3L},
#'                 offering sufficient precision. The maximum limit is twelve
#'                 due to double precision.
#' @param additional_parameters A \code{list} containing additional alternatives for the smoothing spline fit.
#'
#' @description
#' Predict and estimate prediction intervals based on new observations using
#' smoothing splines grouped by \code{comparison}.
#'
#' @details
#' See documentation of \code{predict_smoothing_spline} for detailed information
#' on each of the function arguments.
#'
#'
#' @return A \code{data.table} object comprising approximated prediction interval data based on user inputs.
#' @export
#'
#' @examples print(1)

predict_smoothing_splines <- function(data,
                                      new_data = NULL,
                                      weighted = FALSE,
                                      df = NULL,
                                      lambda = NULL,
                                      df_max = 7.5,
                                      R_ratio = 1,
                                      level = 0.99,
                                      simultaneous = FALSE,
                                      negative_ok = TRUE,
                                      attempt_fast = FALSE,
                                      include_prediction_variance = FALSE,
                                      rounding = 3L,
                                      additional_parameters = list("method" = "loocv",
                                                                   "smudge" = 1.4,
                                                                   "c_star_min" = 0.38,
                                                                   "c_star_max" = 1.50,
                                                                   "tol" = 0.5,
                                                                   "window" = 10,
                                                                   "iter" = 1)){

  # Split input data into a list*
  # WEAKNESS POINT(S):
  # 1. This assumes that data is already a data.table and contains comparison
  data_grouped_by_comparison <- split(data, by = "comparison", keep.by = TRUE, sorted = FALSE)

  # Defaults
  new_data_grouped_by_comparison <- NULL

  # If new_data is NULL, generate new_data based on comparison-wise ranges in data
  # WEAKNESS POINT(S):
  # 2. This assumes that MP_B is part of data
  if(is.null(new_data)){
    new_data_grouped_by_comparison <- lapply(data_grouped_by_comparison, FUN = function(x){
      data.table("comparison" = rep(x$comparison[1], 2e2L),
                 "MP_B" = seq(from = min(x$MP_B, na.rm = TRUE),
                              to = max(x$MP_B, na.rm = TRUE),
                              length.out = 2e2L))
    })
  }
  # If new_data is not NULL, split input new_data into a list (MAKE MORE ROBUST LATER)
  # WEAKNESS POINT(S):
  # 3. This assumes that new_data is already a data.table and contains comparison
  else{
    new_data_grouped_by_comparison <- split(new_data, by = "comparison", keep.by = TRUE, sorted = FALSE)
  }

  # Check if splitted lists have same lengths (equivalent to same numb. of comparisons)
  length_match <- length(data_grouped_by_comparison) == length(new_data_grouped_by_comparison)
  name_match <- all(names(data_grouped_by_comparison) == names(new_data_grouped_by_comparison), na.rm = TRUE)

  if(!length_match){
    stop("data and new_data do not have the same lenghts.")
  }

  else if(!name_match){
    stop("data and new_data do not have the same names or, order of names.")
  }

  # Get comparison names and get corresponding integers
  output_comparison_names <- names(new_data_grouped_by_comparison)
  output_comparison_ints <- 1L:length(new_data_grouped_by_comparison)

  # Checks df argument
  if(!is.null(df)){
    if(length(unique(df)) == 1){
      df <- rep(df[1], length(data_grouped_by_comparison))
    }
    else if(!(length(df) == length(data_grouped_by_comparison))){
      stop("df must be of same length as the number of comparisons in data and new_data")
    }
  }

  # Checks lambda argument
  if(!is.null(lambda)){
    if(length(unique(lambda)) == 1){
      lambda <- rep(lambda[1], length(data_grouped_by_comparison))
    }
    else if(!(length(lambda) == length(data_grouped_by_comparison))){
      stop("lambda must be of same length as the number of comparisons in data and new_data")
    }
  }

  # Run predict_smoothing_spline comparison-wise
  output <- sapply(X = 1:length(data_grouped_by_comparison), FUN = function(x){
    predict_smoothing_spline(data = data_grouped_by_comparison[[x]],
                             new_data = new_data_grouped_by_comparison[[x]],
                             weighted = weighted,
                             df = if(!is.null(df)){df[x]}else{NULL},
                             lambda = if(!is.null(lambda)){lambda[x]}else{NULL},
                             df_max = df_max,
                             R_ratio = R_ratio,
                             level = level,
                             simultaneous = simultaneous,
                             negative_ok = negative_ok,
                             attempt_fast = attempt_fast,
                             include_prediction_variance = include_prediction_variance,
                             rounding = rounding,
                             additional_parameters = additional_parameters)
  }, simplify = FALSE)

  # Convert to data.table object...
  output <- rbindlist(output, idcol = "comparison")
  output$comparison <- output_comparison_names[match(output$comparison, output_comparison_ints)]

  # Return output
  return(output)

}

