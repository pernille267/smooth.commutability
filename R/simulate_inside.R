#' Internal Function for Inside Calculation Given That df_max is Given
#'
#' This function checks if external quality assessment materials are inside estimated prediction intervals given \code{df_max}
#' @keywords internal
simulate_inside_df_max <- function(simulated_cs_data, simulated_eq_data, df_max, level = 0.95){
  # Set default df_0 to NULL
  df_0 <- NULL
  # Calculate interior knots
  interior_knots <- calculate_interior_knots((sort(simulated_cs_data$MP_B) - min(simulated_cs_data$MP_B)) / diff(range(simulated_cs_data$MP_B)))
  # Attempt to fit a smoothing spline using smooth.spline
  ss_fit <- tryCatch(expr = {smooth.spline(x = simulated_cs_data$MP_B, y = simulated_cs_data$MP_A, all.knots = c(0, interior_knots, 1), cv = TRUE, control.spar = list(low = 0.38, high = 1.5))},
                             error = function(e){return(NULL)})
  # If the smooth.spline fit fails, assign lambda to NULL
  if(is.null(ss_fit)){
    lambda_instead <- NULL
    rm(ss_fit)
  }
  # If the smooth.spline fit succeeds, use the lambda
  else{
    lambda_instead <- ss_fit$lambda
    df_0 <- ss_fit$df
    rm(ss_fit)
  }

  # If the df_0 is too large try fitting smooth.spline again using df = df_max
  if((!is.null(df_0)) & df_0 > df_max){
    lambda_instead <- tryCatch(expr = {smooth.spline(x = simulated_cs_data$MP_B, y = simulated_cs_data$MP_A, df = df_max, all.knots = c(0, interior_knots, 1), cv = TRUE)$lambda},
                               error = function(e){return(NULL)})
  }
  else if((!is.null(df_0)) & df_0 <= df_max & df_0 >= 2){
    lambda_instead <- tryCatch(expr = {smooth.spline(x = simulated_cs_data$MP_B, y = simulated_cs_data$MP_A, df = df_0, all.knots = c(0, interior_knots, 1), cv = TRUE)$lambda},
                               error = function(e){return(NULL)})
  }
  else if((!is.null(df_0)) & df_0 < 2){
    lambda_instead <- NULL
  }
  if(is.null(lambda_instead)){
    tryCatch(expr = {
      predict_smoothing_spline(data = simulated_cs_data,
                               new_data = simulated_eq_data,
                               df = NULL,
                               lambda = NULL,
                               df_max = df_max,
                               R_ratio = 1,
                               level = level,
                               rounding = 6L)$inside
    },
    error = function(e){
      return(NA_integer_)
    })
  }
  else{
    tryCatch(expr = {
      predict_smoothing_spline(data = simulated_cs_data,
                               new_data = simulated_eq_data,
                               df = NULL,
                               lambda = lambda_instead,
                               df_max = NULL,
                               R_ratio = 1,
                               level = level,
                               rounding = 6L)$inside
    },
    error = function(e){
      return(NA_integer_)
    })
  }
}

#' Internal Function for Inside Calculation Given That df is Given
#'
#' This function checks if external quality assessment materials are inside estimated prediction intervals given \code{df}
#' @keywords internal
simulate_inside_df <- function(simulated_cs_data, simulated_eq_data, df, level = 0.95){
  interior_knots <- calculate_interior_knots((sort(simulated_cs_data$MP_B) - min(simulated_cs_data$MP_B)) / diff(range(simulated_cs_data$MP_B)))
  lambda_instead <- tryCatch(expr = {smooth.spline(x = simulated_cs_data$MP_B, y = simulated_cs_data$MP_A, df = df, all.knots = c(0, interior_knots, 1), cv = TRUE)$lambda},
                             error = function(e){return(NULL)})
  if(is.null(lambda_instead)){
    tryCatch(expr = {
      predict_smoothing_spline(data = simulated_cs_data,
                               new_data = simulated_eq_data,
                               df = df,
                               lambda = NULL,
                               df_max = NULL,
                               R_ratio = 1,
                               level = level,
                               rounding = 6L)$inside
    },
    error = function(e){
      return(NA_integer_)
    })
  }
  else{
    tryCatch(expr = {
      predict_smoothing_spline(data = simulated_cs_data,
                               new_data = simulated_eq_data,
                               df = NULL,
                               lambda = lambda_instead,
                               df_max = NULL,
                               R_ratio = 1,
                               level = level,
                               rounding = 6L)$inside
    },
    error = function(e){
      return(NA_integer_)
    })
  }
}

#' Simulate Clinical Sample Data and External Quality Assessment Data and Check Whether Inside PI
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#'                   Otherwise, if \code{type} is not included, \code{simulate_eqa_data()} from the \code{fasteqa} package will be triggered, with all relevant
#'                   parameters accepted.
#' @param m A \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#'
#' @details It is important that \code{parameters} holds exactly one value for each simulation parameter. Thus, if \code{parameters} is a \code{data.table} or a \code{data.frame},
#'          it should only have one row. Note that this function only gives one simulation replication. To repeat the simulations multiple times, you can either use the native \code{replicate()}
#'          function, or something similar. Alternatively, you can use the \code{simulate_confidence_level()} function of this package to do the same.
#'
#'
#' @return A \code{m}-dimensional \code{integer} vector with ones and zeros. \code{1} signify that an external quality assessment material is inside the estimated prediction interval, whereas \code{0} signify that it is outside.
#' @export
#'
#' @examples simulate_inside(list(n = 20, R = 3, cvx = 0.04, cvy = 0.01, df = 5))
simulate_inside <- function(parameters, m = 1, level = 0.95){
  SampleID <- MP_B <- NULL
  require_df_max <- FALSE
  if(!any("df" == names(parameters))){
    require_df_max <- TRUE
    df <- NULL
  }
  else if(any("df" == names(parameters))){
    df <- parameters$df
    df_max <- NULL
  }
  if((!any("df_max" == names(parameters))) & require_df_max){
    stop("df_max was not found among the parameters. If df is not part of the parameters, df_max is required.")
  }
  else if((any("df_max" == names(parameters))) & require_df_max){
    df_max <- parameters$df_max
  }

  if(!any("n" == names(parameters))){
    stop("n is missing. It must be part of parameters!")
  }

  if(!is.numeric(m)){
    stop("m is not numeric. m should be an integer larger than or equal to 1.")
  }
  else if(is.numeric(m)){
    if(abs(round(m) - m) > 1e-6){
      m <- round(m)
      warning("m is not an integer. It will be rounded to its nearest integer.")
    }
  }
  if(!is.numeric(level)){
    stop("level is not numeric. level should be a numeric value in the continuous real interval (0, 1).")
  }
  else if(is.numeric(m)){
    if(level <= 0 | level >= 1){
      stop("level should not exceed 1 or be below 0.")
    }
  }

  parameters$n <- parameters$n + m
  if(any("type" == names(parameters))){
    type <- parameters$type
    if(!is.numeric(type)){
      stop("type in parameters is not numeric. It must be an integer and one of the following values: 1, 2, 3.")
    }
    else if(is.numeric(type)){
      if(abs(round(type) - type) > 1e-6){
        type <- round(type)
        if(type > 3L){
          type <- 3
        }
        else if(type < 1){
          type <- 1
        }
        warning("type is not an integer. It will be rounded to its nearest valid integer.")
      }
    }
    if(!any("cve" == names(parameters))){
      parameters$cve <- 0
    }
    simulated_cs_data <- simulate_eqa_data2(parameters = parameters, type = type, AR = FALSE) |> setDT()
    simulated_eq_data <- simulated_cs_data[SampleID %in% sample(x = SampleID[MP_B > min(MP_B) & MP_B < max(MP_B)], size = m, replace = FALSE),]
    simulated_cs_data <- simulated_cs_data[!SampleID %in% simulated_eq_data$SampleID]
    if(sum(is.na(simulated_cs_data$MP_B)) + sum(is.na(simulated_cs_data$MP_A)) > 0){
      simulated_cs_data <- simulate_eqa_data2(parameters = parameters, type = type, AR = FALSE) |> setDT()
      simulated_eq_data <- simulated_cs_data[SampleID %in% sample(x = SampleID[MP_B > min(MP_B) & MP_B < max(MP_B)], size = m, replace = FALSE),]
      simulated_cs_data <- simulated_cs_data[!SampleID %in% simulated_eq_data$SampleID]
      if(sum(is.na(simulated_cs_data$MP_B)) + sum(is.na(simulated_cs_data$MP_A)) > 0){
        return(NA_integer_)
      }
    }

  }
  else{
    simulated_cs_data <- simulate_eqa_data(parameters = parameters)
    simulated_cs_data <- fun_of_replicates(data = simulated_cs_data, fun = "mean") |> setDT()
    simulated_eq_data <- simulated_cs_data[SampleID %in% sample(x = SampleID[MP_B > min(MP_B) & MP_B < max(MP_B)], size = m, replace = FALSE),]
    simulated_cs_data <- simulated_cs_data[!SampleID %in% simulated_eq_data$SampleID]
    if(sum(is.na(simulated_cs_data$MP_B)) + sum(is.na(simulated_cs_data$MP_A)) > 0){
      simulated_cs_data <- simulate_eqa_data(parameters = parameters)
      simulated_cs_data <- fun_of_replicates(data = simulated_cs_data, fun = "mean") |> setDT()
      simulated_eq_data <- simulated_cs_data[SampleID %in% sample(x = SampleID[MP_B > min(MP_B) & MP_B < max(MP_B)], size = m, replace = FALSE),]
      simulated_cs_data <- simulated_cs_data[!SampleID %in% simulated_eq_data$SampleID]
      if(sum(is.na(simulated_cs_data$MP_B)) + sum(is.na(simulated_cs_data$MP_A)) > 0){
        return(NA_integer_)
      }
    }
  }

  if(isFALSE(require_df_max)){
    return(suppressWarnings(simulate_inside_df(simulated_cs_data, simulated_eq_data, df, level)))
  }
  else if(isTRUE(require_df_max)){
    return(suppressWarnings(simulate_inside_df_max(simulated_cs_data, simulated_eq_data, df_max, level)))
  }
  else{
    return(NA_integer_)
  }
}

#' Simulate Clinical Sample Data and External Quality Assessment Data and Check Whether Inside PIs
#'
#' @param N An \code{integer} that represents the number of replicated inside checks.
#' @param m An \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param x Typically an \code{integer}, but not necessarily. Represents the identifier name for the \code{N} inside simulations.
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#'                   Otherwise, if \code{type} is not included, \code{simulate_eqa_data()} from the \code{fasteqa} package will be triggered, with all its choices
#'                   of parameters accepted.
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#' @param attach A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, the relevant simulation parameters are attached to the output for reference.
#' @param replication_id A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, inside check replications for each \code{1:N} is given an identifier in the output. Particularly relevant if \code{m > 1}.
#'
#' @return A \code{data.table} containing the inside replications for the inputs.
#' @export
#'
#' @examples simulate_insides(parameters = list(n = 25, R = 3, cvx = 0.01,
#'  cvy = 0.04, qran = 0.3, qpos = 1, mmax = 5, df = 5))
simulate_insides <- function(N = 1e2L, m = 1, x = 1, parameters, level = 0.95, attach = TRUE, replication_id = TRUE){
  out <- replicate(n = N, expr = simulate_inside(parameters = parameters, m = m, level = level), simplify = FALSE)
  if(m <= 1){
    if(replication_id){
      out <- lapply(out, function(i) data.table(replication = x, inside = i))
    }
    else{
      out <- lapply(out, function(i) data.table(inside = i))
    }

  }
  else{
    if(replication_id){
      out <- sapply(1:length(out), function(i) data.table(replication = x, id = rep(i, length(out[[i]])), SampleID = 1:length(out[[i]]), inside = out[[i]]), simplify = FALSE)
    }
    else{
      out <- sapply(1:length(out), function(i) data.table(SampleID = 1:length(out[[i]]), inside = out[[i]]), simplify = FALSE)
    }
  }
  out <- rbindlist(out)
  if(attach){
    parameters_extended <- lapply(parameters, function(x) rep(x, nrow(out))) |> setDT()
    out <- cbind(parameters_extended, out)
  }
  else{
    return(out)
  }
  return(out)
}

