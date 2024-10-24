#' Internal Function for shifting roles of MP_A and MP_B
#'
#' This function is used internally within \code{fix_one_predictor_and_response()}.
#' @keywords internal
shift_predictor_and_response <- function(data){
  new_data <- copy(data)
  new_data$MP_A <- data$MP_B
  new_data$MP_B <- data$MP_A
  return(new_data)
}

#' Internal Function for shifting roles of MP_A and MP_B if condition is met
#'
#' This function is used internally within \code{fix_predictor_and_response()}.
#' @keywords internal
fix_one_predictor_and_response <- function(data, threshold = 0.5, single_shift = TRUE, force_shift = FALSE){
  if(force_shift){
    if(any("shift" == names(data))){
      if(data$shift[1]){
        if(single_shift){
          return(shift_predictor_and_response(data))
        }
        else{
          return(cbind(shift_predictor_and_response(data), shift = TRUE))
        }
      }
      if(single_shift){
        return(data)
      }
      else{
        return(cbind(data, shift = FALSE))
      }
    }
  }
  impr <- global_precision_estimates(data = data)
  if(impr$lambda < threshold){
    if(single_shift){
      return(shift_predictor_and_response(data))
    }
    else{
      return(cbind(shift_predictor_and_response(data), shift = TRUE))
    }

  }
  if(single_shift){
    return(data)
  }
  else{
    return(cbind(data, shift = FALSE))
  }
}


#' Change roles of MP_A and MP_B based on repeatability variance ratios
#'
#' @param data A \code{list} or \code{data.table} containing columns \code{Comparison},
#'             \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#' @param new_data A \code{list} or \code{data.table} containing columns \code{Comparison},
#'             \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#' @param show_fixes A \code{logical}. If \code{TRUE}, an additonal column is attached to the output.
#'                   This column consists of \code{logical} values. If \code{TRUE}, for a particular
#'                   \code{comparison} element, the roles of MP_A and MP_B are switched for this IVD-MD comparison.
#'
#' @description
#' The smoothing spline model does not account for measurement error in the predictor.
#' To reduce endogeneity-related issues, switching roles of MP_A and MP_B may be a
#' possible solution.
#'
#' @description
#' For each comparison, the roles of MP_A and MP_B are switched if the ratio of
#' repeatability variances of MP_A and MP_B is smaller than 0.5.
#'
#' @return A \code{data.table} containing columns \code{Comparison},
#'         \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#'         If \code{show_fixes} is \code{TRUE}, the output will also contain
#'         the column \code{shift}.
#' @export
#'
#' @examples
#' print(1)
fix_predictor_and_response <- function(data, new_data = NULL, show_fixes = FALSE){
  comparison <- NULL
  if(!isTRUE(is.null(new_data))){
    if(!all(unique(new_data$comparison) == unique(data$comparison))){
      stop("new_data does not contain the same comparisons as data, or is a different order.")
    }
    data <- data[, fix_one_predictor_and_response(.SD, threshold = 0.5, single_shift = FALSE), by = comparison]
    new_data$shift <- data$shift[match(new_data$comparison, data$comparison)]
    new_data <- new_data[, fix_one_predictor_and_response(.SD, threshold = 0.5, single_shift = TRUE, force_shift = TRUE), by = comparison]
    data$comparison <- ifelse(data$shift, sub("(.*) - (.*)", "\\2 - \\1", data$comparison), data$comparison)
    new_data$comparison <- ifelse(new_data$shift, sub("(.*) - (.*)", "\\2 - \\1", new_data$comparison), new_data$comparison)
    if(show_fixes){
      return(list(data = data, new_data = new_data))
    }
    data$shift <- NULL
    new_data$shift <- NULL
    return(list(data = data, new_data = new_data))
  }
  data <- data[, fix_one_predictor_and_response(.SD, threshold = 0.5, single_shift = FALSE), by = comparison]
  data$comparison <- ifelse(data$shift, sub("(.*) - (.*)", "\\2 - \\1", data$comparison), data$comparison)
  if(show_fixes){
    return(data)
  }
  data$shift <- NULL
  return(data)
}
