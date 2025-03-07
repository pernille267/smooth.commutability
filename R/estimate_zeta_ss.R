#' Internal Function for validating entries
#'
#' This function is used internally validate entries in estimate_zeta_ss.
#' @keywords internal
valdiate_estimate_zeta_ss <- function(data, df, weighted, error = FALSE){

  valid <- TRUE
  required_columns_data <- c("SampleID", "ReplicateID", "MP_B", "MP_A")

  # Check if class is correct
  if(!is.data.table(data)){
    if(!is.list(data) && !is.data.frame(data)){
      if(error){
        stop("The data argument is not a data.table, data.frame or list, but a: ",
             class(data))
      }
      valid <- FALSE
    }
  }
  # Check if have the correct names
  missing_columns_data <- required_columns_data[which(!(required_columns_data %in% names(data)))]
  if(length(missing_columns_data) >= 1){
    if(error){
      stop("The data argument does not have the necessary columns. ",
           paste(missing_columns_data, collapse = ", "),
           " are missing.")
    }
    valid <- FALSE
  }

  # Check if df argument is valid
  if(!is.null(df) && (!is.numeric(df) || length(df) != 1)){
    if(error){
      stop("df is expected to be a non-missing numeric value or NULL, but is a ",
           class(df),
           ".")
    }
    valid <- FALSE
  }
  if(!is.null(df) && is.numeric(df) && df >= length(data$SampleID)){
    if(error){
      stop("df cannot exceed the number of unique samples in data.")
    }
  }

  # Check if weighted argument is valid
  if(is.null(weighted)){
    if(error){
      stop("weighted is expected to be a non-missing logical value, ",
           "but is NULL.")
    }
    valid <- FALSE
  }
  else if(!is.logical(weighted)){
    if(error){
      stop("weighted is expected to be a non-missing logical value, ",
           "but is a ",
           class(weighted))
    }
    valid <- FALSE
  }
  else if(length(weighted) != 1){
    if(error){
      stop("weighted is expected to be a non-missing logical value, ",
           "but is a logical vector of length.",
           length(weighted))
    }
    valid <- FALSE
  }
  else if(is.na(weighted)){
    if(error){
      stop("weighted is expected to be a non-missing logical value, ",
           "but is a NA value.")
    }
    valid <- FALSE
  }
  return(valid)
}


#' Estimate \eqn{\zeta} Using Smoothing Splines
#'
#' @param data A \code{list} or \code{data.table} object. Must contain columns,
#'             \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#' @param df A \code{double} between \code{2} and \code{n}, where \code{n} is the
#'           number of unique elements of \code{SampleID}. See details.
#' @param weighted A \code{logical} value. If \code{TRUE}, iteratively estimated
#'                 weights are used to fit the smoothing spline.
#' @param mor A \code{logical} value.
#' @param na_rm A \code{logical} value. Should NA values be removed before
#'              estimating \eqn{\zeta}?
#'
#' @details
#' Estimates \eqn{\zeta} using the plug-in estimator
#'
#' \eqn{\hat{\zeta} = \frac{\hat{\sigma}^2 \cdot \frac{n + \mathrm{df}(c)}{n}}{\hat{\sigma}_v^2 + \hat{\sigma}_h^2 \frac{1}{N} \sum_{i=1}^{n}\sum_{r=1}^{R_i}[\hat{f}(x_{ir})]^2}}
#'
#'
#' @return
#' A \code{list} with one variable named \code{zeta}. This value is a \code{double}
#' that signifies the smoothing spline estimate of \eqn{\zeta} based on \code{df}.
#' @export
#'
#' @examples print(1)

estimate_zeta_ss <- function(data, df = NULL, weighted = FALSE, mor = FALSE, na_rm = TRUE){

  # Initialization
  R <- 1

  # Validate entries
  valid_entries <- valdiate_estimate_zeta_ss(data = data,
                                             df = df,
                                             weighted = weighted,
                                             error = TRUE)

  # If not valid entries => return NA
  if(!valid_entries){
    return(list(zeta = NA_real_))
  }

  # Given that valid_entries must be TRUE here, we can convert data to data.table
  data <- as.data.table(data)

  # If df is 2, we use estimate_zeta() method
  if(!is.null(df) && df <= 2){
    if(!mor){
      return(estimate_zeta_ols(data = data))
    }
  }

  # Glocal precision estimates
  impr <- global_precision_estimates(data = data)

  # Check global precision estimates
  # global_precision_estimates may result in unreliable estimates
  # if CV_A -> 0 or CV_B -> 0. These next lines of code are supposed to
  # guard against these cases.
  if(impr$Var_B <= 0 + .Machine$double.eps){
    impr$Var_B <- 0
    impr$CV_B <- 0
    impr$lambda <- .Machine$double.xmax
  }
  else if(impr$lambda <= 0 + .Machine$double.eps){
    impr$Var_B <- 0
    impr$CV_B <- 0
    impr$lambda <- .Machine$double.eps
  }

  # If na_rm = TRUE, remove present NA values
  if(na_rm){
    not_na_ids <- (!is.na(data$MP_B)) & (!is.na(data$MP_A))
    data <- data[not_na_ids, ]
  }

  # Shift roles of predictor and response if lambda < 0.5
  axis_variables <- c("MP_B", "MP_A", "Var_B", "Var_A")
  if(impr$lambda < 0.5){
    axis_variables <- c("MP_A", "MP_B", "Var_A", "Var_B")
  }

  if(mor){
    R <- count_samplewise_replicates(data = data,
                                     summary = "mean",
                                     invalid_NA = TRUE,
                                     silence = 1)$R_i
    if(is.null(R)){
      stop("R was not correctly calculated. Check count_samplewise_replicates().")
    }
    data <- data[, fun_of_replicates(.SD)]
  }

  # Order according to predictor values
  setorderv(data, axis_variables[1], order = 1)

  # Ordered (x_i, y_i, w_i)
  x <- data[[axis_variables[1]]]
  y <- data[[axis_variables[2]]]
  w <- 1

  # Extract repeatability variances
  var_h <- impr[[axis_variables[3]]] / R
  var_v <- impr[[axis_variables[4]]] / R

  # Get weights
  if(weighted){
    w <- "estimate"
  }

  # Fit the smoothing spline
  ss_fit <- smoothing_spline(data = data,
                             weights = w,
                             df = df,
                             df_max = NULL,
                             attempt_fast = FALSE,
                             na_rm = FALSE)#,
                             #additional_parameters = ...)

  numerator <- sum(ss_fit$var_eps / ss_fit$weights / ss_fit$n)  * (ss_fit$n + ss_fit$df) / ss_fit$n
  denominator <- var_v + var_h * sum(ss_fit$derivatives^2 * ss_fit$weights / sum(ss_fit$weights))
  zeta <- numerator / denominator
  return(list("zeta" = zeta))
}

#' Estimate Comparison-wise \eqn{\zeta} using Smoothing Splines
#'
#' @param data A \code{data.table} object. Must contain columns \code{comparison},
#'             \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}
#' @param df A \code{double} between \code{2} and \code{n}, where \code{n} is the
#'           number of unique elements of \code{SampleID}. Otherwise, it must be a
#'           \code{numeric} vector of length equal to the number of unique elements
#'           in \code{comparison}.
#' @param weighted A \code{logical} value. If \code{TRUE}, iteratively estimated
#'                 weights are used to fit the smoothing spline.
#' @param na_rm A \code{logical} value. Should NA values be removed before
#'              estimating \eqn{\zeta}?
#'
#' @return
#' A \code{data.table} with two variables. The first variable is the \code{comparison}
#' variable. The second variable is the corresponding estimate of \eqn{\zeta},
#' based on \code{df} and \code{weighted}.
#' @export
#'
#' @examples print(1)

estimate_zetas_ss <- function(data, df = NULL, weighted = FALSE, na_rm = TRUE){

  # Initialization
  required_columns_data <- c("comparison",
                             "SampleID",
                             "ReplicateID",
                             "MP_B",
                             "MP_A")
  valid_df <- TRUE
  valid_weighted <- TRUE
  must_split <- TRUE

  # Checking validity of data
  if(!is.data.table(data) & !is.data.frame(data) & !is.list(data)){
    stop("data is not a data.table, data.frame or list, but a: ",
         class(data)[1])
  }
  else if(!is.data.table(data)){
    if(is.data.frame(data)){
      data <- as.data.table(data)
    }
    else if(is.list(data)){
      if(length(required_columns_data[which(!(required_columns_data %in% names(data)))]) == 0){
        data <- as.data.table(data)
      }
      else if(length(required_columns_data[which(!(required_columns_data %in% names(data)))]) == 5){
        if(is.data.table(data[[1]])){
          if(length(required_columns_data[which(!(required_columns_data %in% names(data[[1]])))]) == 4){
            must_split <- FALSE
          }
          else{
            stop("data seems to be already a list grouped by comparison, but ",
                 "the first element does not have the correct column names.")
          }
        }
        else{
          stop("data seems to be already a list grouped by comparison, but ",
               "the first element is not a data.table.")
        }
      }
    }
  }
  if(!is.data.table(data)){
    stop("data is attempted to be converted to data.table, but is still",
         " not a data.table, but a ",
         class(data)[1])
  }
  missing_columns_data <- required_columns_data[which(!(required_columns_data %in% names(data)))]
  if(length(missing_columns_data) >= 1){
    stop("The data argument does not have the necessary columns. ",
         paste(missing_columns_data, collapse = ", "),
         " are missing.")
  }

  # Split data by comparison column
  if(must_split){
    data_grouped_by_comparison <- split(copy(data),
                                        by = "comparison",
                                        keep.by = TRUE)
  }
  else{
    data_grouped_by_comparison <- copy(data)
  }

  # Get unique comparison names
  output_comparison_names <- names(data_grouped_by_comparison)

  # Checks if some of the names are invalid
  if(any(required_columns_data %in% output_comparison_names)){
    stop("The splitted data contain some invalid names: ",
         paste(required_columns_data[required_columns_data %in% output_comparison_names],
               collapse = ", "),
         ".")
  }

  # Get number of unique comparisons in data
  nc <- length(output_comparison_names)
  # Get integer represetnation of each unique comparison
  output_comparison_ints <- 1L:nc

  # Check validity of df
  if(!is.null(df)){
    if(!is.numeric(df)){
      stop("df must be NULL, a single double value or a numeric vector ",
           "of same length as the number of unique comparisons in data.")
    }
    else if(length(df) == 0){
      stop("df must be NULL, a single double value or a numeric vector ",
           "of same length as the number of unique comparisons in data.")
    }
    else if(length(df) >= 2 & length(df) < nc){
      stop("length of df is larger than 1, but not equal to the number ",
           "of unique comparisons in data.")
    }
    else if(length(df) == 1){
      df <- rep(df, nc)
    }
  }

  # Estimate zeta values for each unique comparison
  zeta_grouped_by_comparison <- sapply(X = 1:nc,
                                       FUN = function(x){
                                         estimate_zeta_ss(data = data_grouped_by_comparison[[x]],
                                                          df = df[x],
                                                          weighted = weighted,
                                                          na_rm = na_rm)
                                       }, simplify = FALSE)

  zeta_grouped_by_comparison <- rbindlist(zeta_grouped_by_comparison,
                                          idcol = "comparison")
  zeta_grouped_by_comparison$comparison <- output_comparison_names[match(zeta_grouped_by_comparison$comparison,
                                                                         output_comparison_ints)]

  return(zeta_grouped_by_comparison)

}

