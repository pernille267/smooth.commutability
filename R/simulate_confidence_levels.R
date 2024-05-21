#' Internal Function for Checking Input Parameters
#'
#' This function is used internally for checking if the inputs are valid
#' @keywords internal
parameter_checks <- function(N, n, m, level, attach, global_mean, all_inside, calculate_sem){
  if(isTRUE(is.na(N)) | !is.numeric(N)){
    stop("N is not numeric or is missing. It should be an integer larger than or equal to 2.")
  }
  else if(isTRUE(is.na(n)) | !is.numeric(n)){
    stop("n is not numeric or is missing. It should be an integer larger than or equal to 1.")
  }
  else if(isTRUE(is.na(m)) | !is.numeric(m)){
    stop("m is not numeric or is missing. It should be an integer larger than or equal to 1.")
  }
  else if(isTRUE(is.na(level)) | !is.numeric(level)){
    stop("level is not numeric or is missing. It should be a double in the real-open interval (0, 1).")
  }
  else if(isTRUE(is.na(attach)) | !is.logical(attach)){
    stop("attach is not a non-missing logical value. It should be either TRUE or FALSE.")
  }
  else if(isTRUE(is.na(global_mean)) | !is.logical(global_mean)){
    stop("global_mean is not a non-missing logical value. It should be either TRUE or FALSE.")
  }
  else if(isTRUE(is.na(all_inside)) | !is.logical(all_inside)){
    stop("all_inside is not a non-missing logical value. It should be either TRUE or FALSE.")
  }
  else if(isTRUE(is.na(calculate_sem)) | !is.logical(calculate_sem)){
    stop("calculate_sem is not a non-missing logical value. It should be either TRUE or FALSE.")
  }
  else{
    N <- N[1]
    n <- n[1]
    m <- m[1]
    level <- level[1]
    attach <- attach[1]
    global_mean <- global_mean[1]
    all_inside <- all_inside[1]
    calculate_sem <- calculate_sem[1]
    if(abs(round(N) - N) > 1e-6){
      stop("N is not an integer")
    }
    else if(N < 2){
      stop("N < 2, but N >= 2 is required.")
    }
    else if(N > 1e6L){
      stop("N > 1e6L, but this will take forever to run. Was this your intention?")
    }

    if(abs(round(n) - n) > 1e-6){
      stop("n is not an integer")
    }
    else if(n < 1){
      stop("n < 1, but n >= 1 is required.")
    }
    else if(n > 1e3L){
      stop("n > 1e3L, but this will take forever to run. Was this your intention?")
    }

    if(m > 1e3L){
      stop("m > 1e3L, but this choice is insane. Why??")
    }

    if(level <= 0 | level >= 1){
      stop("level is outside the real-open interval (0, 1). Why??")
    }
  }
}

#' Simulate empirical confidence level for a set of parameters
#'
#' @param N An \code{integer} larger than or equal to \code{1L} that represents the number of replicated inside checks.
#' @param n An \code{integer} larger than or equal to \code{1L} that represents the number of replicated confidence levels for the particular \code{parameters}.
#' @param m An \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#'                   Otherwise, if \code{type} is not included, \code{simulate_eqa_data()} from the \code{fasteqa} package will be triggered, with all its choices
#'                   of parameters accepted.
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#' @param attach A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, the relevant simulation parameters are attached to the output for reference.
#' @param global_mean A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, empirical confidence level is estimated taking the mean overall inside checks.
#' @param all_inside A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE} and \code{m > 1}, family-wise emprical confidence level is estimated.
#' @param calculate_sem A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, standard error of the empirical confidence level(s) are calculated.
#'
#' @return A \code{data.table} with \code{n} rows, containing the \code{1:n} empirical confidence level estimates for the particular \code{parameters}.
#' @export
#'
#' @examples print(1)
simulate_confidence_level <- function(N = 1e2L, n = 1, m = 1, parameters, level = 0.95, attach = TRUE, global_mean = FALSE, all_inside = FALSE, calculate_sem = TRUE){
  parameter_checks(N, n, m, level, attach, global_mean, all_inside, calculate_sem)
  inside <- NULL
  if(n > 1){
    out <- sapply(X = 1:n, function(i) simulate_insides(N = N, m = m, x = i, parameters = parameters, level = level, attach = attach, replication_id = TRUE), simplify = FALSE)
    out <- rbindlist(out)
    not_group <- c("inside", "SampleID")
    if(m <= 1){
      print("m <= 1")
      all_inside <- FALSE
      global_mean <- FALSE
    }
    if(all_inside){
      print("all_inside")
      global_mean <- FALSE
      out <- out[, list(inside = as.integer(identical(inside, rep(1L, length(inside))))), by = setdiff(names(out), not_group)]
      print(out)
    }

    if(global_mean){
      if(m > 1){
        not_group <- c("inside", "SampleID", "id")
      }
      out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                        number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), not_group)]
      return(out)
    }
    if((!global_mean) & (!all_inside)){
      print("not global_mean & not all_inside")
      if(m > 1){
        out <- out[, list(inside = mean(inside, na.rm = TRUE)), by = setdiff(names(out), not_group)]
      }
    }

    if(calculate_sem){
      if(m > 1){
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside)),
                          sem = sqrt(mean(inside, na.rm = TRUE) * (1 - mean(inside, na.rm = TRUE)) / sum(!is.na(inside)))), by = setdiff(names(out), c("inside", "id"))]
      }
      else{
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside)),
                          sem = sqrt(mean(inside, na.rm = TRUE) * (1 - mean(inside, na.rm = TRUE)) / sum(!is.na(inside)))), by = setdiff(names(out), c("inside"))]
      }

    }
    else{
      if(m > 1){
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), c("inside", "id"))]
      }
      else{
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), c("inside"))]
      }

    }
    return(out)
  }
  else if(n == 1){
    out <- simulate_insides(N = N, m = m, x = 1, parameters = parameters, level = level, attach = attach, replication_id = TRUE)
    out$replication <- NULL
    not_group <- c("inside", "SampleID")
    if(m <= 1){
      print("m <= 1")
      all_inside <- FALSE
      global_mean <- FALSE
    }
    if(all_inside){
      print("all_inside")
      global_mean <- FALSE
      out <- out[, list(inside = as.integer(identical(inside, rep(1L, length(inside))))), by = setdiff(names(out), not_group)]
      print(out)
    }

    if(global_mean){
      if(m > 1){
        not_group <- c("inside", "SampleID", "id")
      }
      out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                        number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), not_group)]
      return(out)
    }
    if((!global_mean) & (!all_inside)){
      print("not global_mean & not all_inside")
      if(m > 1){
        out <- out[, list(inside = mean(inside, na.rm = TRUE)), by = setdiff(names(out), not_group)]
      }
    }

    if(calculate_sem){
      if(m > 1){
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside)),
                          sem = sqrt(mean(inside, na.rm = TRUE) * (1 - mean(inside, na.rm = TRUE)) / sum(!is.na(inside)))), by = setdiff(names(out), c("inside", "id"))]
      }
      else{
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside)),
                          sem = sqrt(mean(inside, na.rm = TRUE) * (1 - mean(inside, na.rm = TRUE)) / sum(!is.na(inside)))), by = setdiff(names(out), c("inside"))]
      }

    }
    else{
      if(m > 1){
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), c("inside", "id"))]
      }
      else{
        out <- out[, list(confidence_level = mean(inside, na.rm = TRUE),
                          number_of_observations = sum(!is.na(inside))), by = setdiff(names(out), c("inside"))]
      }

    }
    return(out)
  }
  else{
    stop("NO!")
  }
}

#' Simulate empirical confidence levels for several sets of parameters
#'
#' @param N An \code{integer} larger than or equal to \code{1L} that represents the number of replicated inside checks.
#' @param n An \code{integer} larger than or equal to \code{1L} that represents the number of replicated confidence levels for the particular \code{parameters}.
#' @param m An \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param parameters A \code{data.table} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#'                   Otherwise, if \code{type} is not included, \code{simulate_eqa_data()} from the \code{fasteqa} package will be triggered, with all its choices
#'                   of parameters accepted.
#' @param level A \code{double} that is larger than 0, but smaller than 1. Represents the confidence level of the estimated prediction intervals for the smoothing spline fit.
#' @param attach A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, the relevant simulation parameters are attached to the output for reference.
#' @param global_mean A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, empirical confidence level is estimated taking the mean overall inside checks.
#' @param all_inside A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE} and \code{m > 1}, family-wise emprical confidence level is estimated.
#' @param calculate_sem A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, standard error of the empirical confidence level(s) are calculated.
#' @param percent A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, all numbers that are possible to present as percentages are presented as percentages.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#'
#' @return A \code{data.table} with \code{n} rows, containing the \code{1:n} empirical confidence level estimates for the particular \code{parameters}.
#' @export
#'
#' @examples print(1)
simulate_confidence_levels <- function(N = 1e2L, n = 1, m = 1, parameters, level = 0.95, attach = TRUE, global_mean = FALSE, all_inside = FALSE, calculate_sem = TRUE, percent = FALSE, parallel = TRUE){
  parameters_list <- split(x = parameters, f = 1:nrow(parameters))

  if(isTRUE(parallel)){
    new_enviroment <- new.env()
    new_enviroment$N <- N
    new_enviroment$n <- n
    new_enviroment$m <- m
    new_enviroment$parameters_list <- parameters_list
    new_enviroment$level <- level
    new_enviroment$attach <- attach
    new_enviroment$global_mean <- global_mean
    new_enviroment$all_inside <- all_inside
    new_enviroment$calculate_sem <- calculate_sem
    new_enviroment$simulate_confidence_level <- simulate_confidence_level
    new_enviroment$simulate_insides <- simulate_insides
    new_enviroment$simulate_inside <- simulate_inside
    new_enviroment$simulate_inside_df <- simulate_inside_df
    new_enviroment$simulate_inside_df_max <- simulate_inside_df_max
    new_enviroment$parameter_checks <- parameter_checks

    cl <- makeCluster(detectCores() - 1)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl = cl, expr = {library(data.table);library(fasteqa)})
    clusterExport(cl = cl, varlist = c("N", "n", "m", "parameters_list", "level", "attach", "global_mean", "all_inside", "calculate_sem", "simulate_confidence_level", "simulate_insides", "simulate_inside", "simulate_inside_df", "simulate_inside_df_max", "parameter_checks"), envir = new_enviroment)
    out <- pblapply(X = parameters_list, FUN = function(x) simulate_confidence_level(N = N, n = n, m = m, parameters = x, level = level, attach = attach, global_mean = global_mean, all_inside = all_inside, calculate_sem = calculate_sem), cl = cl)
    #stopCluster(cl = cl)
    out <- rbindlist(out)
  }
  else if(!isTRUE(parallel)){
    out <- pblapply(X = parameters_list, FUN = function(x) simulate_confidence_level(N = N, n = n, m = m, parameters = x, level = level, attach = attach, global_mean = global_mean, all_inside = all_inside, calculate_sem = calculate_sem))
    out <- rbindlist(out)
  }

  if(isTRUE(percent)){
    if(any("cvx" == names(out))){
      out$cvx <- out$cvx * 100
    }
    if(any("cvy" == names(out))){
      out$cvy <- out$cvy * 100
    }
    if(any("cve" == names(out))){
      out$cve <- out$cve * 100
    }
    if(any("confidence_level" == names(out))){
      out$confidence_level <- out$confidence_level * 100
    }
    if(any("sem" == names(out))){
      out$sem <- out$sem * 100
    }
  }

  return(out)
}

