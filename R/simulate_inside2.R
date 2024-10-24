#' Simulate Clinical Sample Data and External Quality Assessment Data and Check Whether Inside PI
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param m A \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#' @param extrapolate_ok A \code{logical} value. If \code{TRUE}, test data can be simulated beyond the empirical range of the clinical sample data.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#'
#' @return An \code{integer} vector of length \code{m}.
#' @export
#'
#' @examples print(1)
simulate_inside2 <- function(parameters, m = 1, level = 0.95, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE){
  type <- NULL
  shift_roles <- shift & all(c("cvx", "cvy") %in% names(parameters)) & ((parameters$cvy / parameters$cvx) ** 2) < 0.5
  if("type" %in% names(parameters)){
    type <- parameters$type[1]
  }
  else{
    type <- 0
  }
  simulated_cs_data <- simulate_eqa_data2(parameters = parameters, type = type, AR = FALSE, include_parameters = TRUE, shift = shift_roles)
  eq_parameters <- simulated_cs_data$parameters
  eq_parameters$n <- m
  simulated_cs_data <- simulated_cs_data$simulated_data
  if(!extrapolate_ok){
    eq_parameters$cil <- min(simulated_cs_data$MP_B, na.rm = TRUE)
    eq_parameters$ciu <- max(simulated_cs_data$MP_B, na.rm = TRUE)
  }
  simulated_eq_data <- simulate_eqa_data2(parameters = eq_parameters, type = type, AR = FALSE, include_parameters = FALSE, shift = shift_roles)
  output <- tryCatch(expr = predict_smoothing_spline(data = simulated_cs_data,
                                                     new_data = simulated_eq_data,
                                                     weight_data = NULL,
                                                     df = parameters$df,
                                                     lambda = parameters$lambda,
                                                     df_max = parameters$df_max,
                                                     R_ratio = 1,
                                                     level = level,
                                                     negative_ok = TRUE,
                                                     attempt_fast = attempt_fast,
                                                     rounding = 3L)$inside,
                     error = function(e) rep(NA_integer_, m))
  return(output)
}

#' Internal Function for Inside Calculation
#'
#' This function checks if external quality assessment materials are inside estimated prediction intervals
#' @keywords internal
simulate_insides2_internal <- function(parameters, N = 100, m = 1, level = 0.95, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE){
  mid <- rep(1:m, times = N)
  nid <- rep(1:N, each = m)
  ins <- replicate(n = N, expr = simulate_inside2(parameters = parameters, m = m, level = level, extrapolate_ok = extrapolate_ok, shift = shift, attempt_fast = attempt_fast), simplify = TRUE)
  spa <- lapply(parameters, rep, times = N * m)
  out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins))
  return(out)
}

#' Simulate Clinical Sample Data and External Quality Assessment Data and Check Whether Inside PIs
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param N An \code{integer} that represents the number of replicated inside checks.
#' @param m An \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#' @param extrapolate_ok A \code{logical} value. If \code{TRUE}, test data can be simulated beyond the empirical range of the clinical sample data.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#'
#' @return A \code{list}.
#' @export
#'
#' @examples print(1)
simulate_insides2 <- function(parameters, N = 100, m = 1, level = 0.95, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE, parallel = TRUE, max_cores = 4){
  multiple_pars <- FALSE
  must_split <- FALSE
  if(is.data.table(parameters)){
    if(nrow(parameters) > 1){
      multiple_pars <- TRUE
      must_split <- TRUE
    }
  }
  else if(is.list(parameters)){
    must_split <- !all(sapply(parameters, is.data.table, USE.NAMES = FALSE))
  }
  if(multiple_pars){
    if(must_split){
      parameters_list <- split(parameters, f = 1:nrow(parameters))
    }
    else{
      parameters_list <- parameters
    }

    if(isTRUE(parallel)){
      new_enviroment <- new.env()
      new_enviroment$N <- N
      new_enviroment$m <- m
      new_enviroment$parameters_list <- parameters_list
      new_enviroment$level <- level
      new_enviroment$extrapolate_ok <- extrapolate_ok
      new_enviroment$shift <- shift
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$simulate_insides2_internal <- simulate_insides2_internal
      new_enviroment$simulate_inside2 <- simulate_inside2
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterEvalQ(cl = cl, expr = {library(data.table)})
      clusterExport(cl = cl, varlist = c("N", "m", "parameters_list", "level", "extrapolate_ok", "shift", "attempt_fast", "simulate_eqa_data2", "simulate_insides2_internal", "simulate_inside2"), envir = new_enviroment)
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_insides2_internal(parameters = x, N = N, m = m, level = level, extrapolate_ok = extrapolate_ok, shift = shift, attempt_fast = attempt_fast), cl = cl)
      out <- rbindlist(out)
    }
    else if(!isTRUE(parallel)){
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_insides2_internal(parameters = x, N = N, m = m, level = level, extrapolate_ok = extrapolate_ok, shift = shift, attempt_fast = attempt_fast))
      out <- rbindlist(out)
    }
  }
  else if(!multiple_pars){
    if(isTRUE(parallel)){
      new_enviroment <- new.env()
      new_enviroment$N <- N
      new_enviroment$m <- m
      new_enviroment$parameters <- parameters
      new_enviroment$level <- level
      new_enviroment$extrapolate_ok <- extrapolate_ok
      new_enviroment$shift <- shift
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$simulate_insides2_internal <- simulate_insides2_internal
      new_enviroment$simulate_inside2 <- simulate_inside2
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterEvalQ(cl = cl, expr = {library(data.table);library(fasteqa)})
      clusterExport(cl = cl, varlist = c("N", "m", "parameters", "level", "extrapolate_ok", "shift", "attempt_fast", "simulate_insides2_internal", "simulate_inside2", "simulate_eqa_data2"), envir = new_enviroment)
      mid <- rep(1:m, times = N)
      nid <- rep(1:N, each = m)
      spa <- lapply(parameters, rep, times = N * m)
      ins <- pbreplicate(n = N, expr = simulate_inside2(parameters = parameters, m = m, level = level, extrapolate_ok = extrapolate_ok, shift = shift, attempt_fast = attempt_fast), cl = cl, simplify = TRUE)
      out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins)) |> setDT()

    }
    else if(!isTRUE(parallel)){
      mid <- rep(1:m, times = N)
      nid <- rep(1:N, each = m)
      spa <- lapply(parameters, rep, times = N * m)
      ins <- pbreplicate(n = N, expr = simulate_inside2(parameters = parameters, m = m, level = level, extrapolate_ok = extrapolate_ok, shift = shift, attempt_fast = attempt_fast), simplify = TRUE)
      out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins)) |> setDT()
    }
  }
  return(out)
}







