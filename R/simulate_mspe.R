#' Simulate Prediction Errors
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param m A \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param extrapolate_ok A \code{logical} value. If \code{TRUE}, test data can be simulated beyond the empirical range of the clinical sample data.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#'
#' @return An \code{numeric} vector of length \code{m}.
#' @export
#'
#' @examples print(1)
simulate_pe <- function(parameters,
                        m = 1,
                        extrapolate_ok = FALSE,
                        shift = FALSE,
                        attempt_fast = FALSE){
  # Initialize
  type <- NULL
  obs_tau <- NULL
  additional_parameters <- list("method" = "gcv",
                                "c_star_min" = 0.38,
                                "c_star_max" = 1.50,
                                "smudge" = 1.2,
                                "tol" = 0.5,
                                "window" = 10,
                                "iter" = 1)


  shift_roles <- shift & all(c("cvx", "cvy") %in% names(parameters)) & ((parameters$cvy / parameters$cvx) ** 2) < 0.5
  if("type" %in% names(parameters)){
    type <- parameters$type[1]
  }
  else{
    type <- 0
  }
  if("obs_tau" %in% names(parameters)){
    obs_tau <- parameters$obs_tau
    parameters$obs_tau <- NULL
  }
  if("method" %in% names(parameters)){
    additional_parameters$method <- parameters$method
  }
  if("c_star_min" %in% names(parameters)){
    additional_parameters$c_star_min <- parameters$c_star_min
  }
  if("c_star_max" %in% names(parameters)){
    additional_parameters$c_star_max <- parameters$c_star_max
  }
  if("smudge" %in% names(parameters)){
    additional_parameters$smudge <- parameters$smudge
  }
  simulated_cs_data <- simulate_eqa_data2(parameters = parameters, type = type, AR = FALSE, include_parameters = TRUE, shift = shift_roles)
  if(!is.null(obs_tau)){
    eq_parameters <- simulated_cs_data$parameters
    eq_parameters$n <- m
    eq_parameters$obs_tau <- obs_tau
    simulated_cs_data <- simulated_cs_data$simulated_data
    simulated_eq_data <- simulate_eqa_data2(parameters = eq_parameters, type = type, AR = FALSE, include_parameters = FALSE, shift = shift_roles)
    output_ss <- tryCatch(expr = predict_smoothing_spline(data = simulated_cs_data,
                                                          new_data = simulated_eq_data,
                                                          weighted = FALSE,
                                                          df = if(any("df" == names(parameters))){parameters$df}else{NULL},
                                                          lambda = parameters$lambda,
                                                          df_max = parameters$df_max,
                                                          R_ratio = 1,
                                                          level = 0.95,
                                                          negative_ok = TRUE,
                                                          attempt_fast = attempt_fast,
                                                          rounding = 3L,
                                                          additional_parameters = additional_parameters)$prediction,
                          error = function(e) rep(NA_real_, m))
    output_ss <- simulated_eq_data$MP_A - output_ss

    return(output_ss)
  }
  eq_parameters <- simulated_cs_data$parameters
  eq_parameters$n <- m
  if(!extrapolate_ok){
    eq_parameters$cil <- min(simulated_cs_data$simulated_data$MP_B, na.rm = TRUE)
    eq_parameters$ciu <- max(simulated_cs_data$simulated_data$MP_B, na.rm = TRUE)
    if(eq_parameters$ciu <= eq_parameters$cil){
      return(rep(NA_integer_, m))
    }
    eq_parameters$cvx <- NULL # Remove because we need to force use of sdx instead
    eq_parameters$cvy <- NULL # Remove because we need to force use of sdy instead
  }
  simulated_cs_data <- simulated_cs_data$simulated_data
  simulated_eq_data <- simulate_eqa_data2(parameters = eq_parameters, type = type, AR = FALSE, include_parameters = FALSE, shift = shift_roles)
  output_ss <- tryCatch(expr = predict_smoothing_spline(data = simulated_cs_data,
                                                        new_data = simulated_eq_data,
                                                        weighted = FALSE,
                                                        df = if(any("df" == names(parameters))){parameters$df}else{NULL},
                                                        lambda = parameters$lambda,
                                                        df_max = parameters$df_max,
                                                        R_ratio = 1,
                                                        level = 0.95,
                                                        negative_ok = TRUE,
                                                        attempt_fast = attempt_fast,
                                                        rounding = 3L,
                                                        additional_parameters = additional_parameters)$prediction,
                        error = function(e) rep(NA_real_, m))
  output_ss <- simulated_eq_data$MP_A - output_ss
  return(output_ss)
}

#' Internal Function for Prediction Error Calculation
#'
#' This function replicates prediction error calculation \code{N} times for \code{parameters}
#' @keywords internal
simulate_pes_internal <- function(parameters,
                                  N = 100,
                                  m = 1,
                                  extrapolate_ok = FALSE,
                                  shift = FALSE,
                                  attempt_fast = FALSE){
  mid <- rep(1:m, times = N)
  nid <- rep(1:N, each = m)
  pes <- replicate(n = N, expr = simulate_pe(parameters = parameters,
                                             m = m,
                                             extrapolate_ok = extrapolate_ok,
                                             shift = shift,
                                             attempt_fast = attempt_fast), simplify = TRUE) |> as.vector()
  spa <- lapply(parameters, rep, times = N * m)
  out <- c(spa, list("m" = mid, "x" = nid, "pe" = pes))
  return(out)
}

#' Replicate Calculation of Prediction Error for Each Parameter Combination
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param N An \code{integer} that represents the number of replicated MSPE calculations
#'          based on each parameter combination.
#' @param m An \code{integer} larger than or equal to \code{1L}. How many external quality assessment materials should be simulated?
#' @param extrapolate_ok A \code{logical} value. If \code{TRUE}, test data can be simulated beyond the empirical range of the clinical sample data.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here. If you do not want reproducibility, set this to \code{NULL}.
#'
#' @return A \code{data.table}.
#' @export
#'
#' @examples print(1)
simulate_pes <- function(parameters, N = 100, m = 1, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE, parallel = TRUE, max_cores = 4, seed = 99){
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
      new_enviroment$extrapolate_ok <- extrapolate_ok
      new_enviroment$shift <- shift
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$simulate_pes_internal <- simulate_pes_internal
      new_enviroment$simulate_pe <- simulate_pe
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "m",
                                         "parameters_list",
                                         "extrapolate_ok",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "simulate_eqa_data2",
                                         "simulate_pes_internal",
                                         "simulate_pe"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_pes_internal(parameters = x,
                                                                                   N = N,
                                                                                   m = m,
                                                                                   extrapolate_ok = extrapolate_ok,
                                                                                   shift = shift,
                                                                                   attempt_fast = attempt_fast), cl = cl)
      out <- rbindlist(out)
    }
    else if(!isTRUE(parallel)){
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_pes_internal(parameters = x,
                                                                                   N = N,
                                                                                   m = m,
                                                                                   extrapolate_ok = extrapolate_ok,
                                                                                   shift = shift,
                                                                                   attempt_fast = attempt_fast))
      out <- rbindlist(out)
    }
  }
  else if(!multiple_pars){
    if(isTRUE(parallel)){
      new_enviroment <- new.env()
      new_enviroment$N <- N
      new_enviroment$m <- m
      new_enviroment$parameters <- parameters
      new_enviroment$extrapolate_ok <- extrapolate_ok
      new_enviroment$shift <- shift
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$simulate_pe <- simulate_pe
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "m",
                                         "parameters",
                                         "extrapolate_ok",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "simulate_pe",
                                         "simulate_eqa_data2"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      mid <- rep(1:m, times = N)
      nid <- rep(1:N, each = m)
      spa <- lapply(parameters, rep, times = N * m)
      pes <- pbreplicate(n = N, expr = simulate_pe(parameters = parameters,
                                                   m = m,
                                                   extrapolate_ok = extrapolate_ok,
                                                   shift = shift,
                                                   attempt_fast = attempt_fast), cl = cl, simplify = TRUE) |> as.vector()
      out <- c(spa, list("m" = mid, "x" = nid, "pe" = pes)) |> setDT()

    }
    else if(!isTRUE(parallel)){
      mid <- rep(1:m, times = N)
      nid <- rep(1:N, each = m)
      spa <- lapply(parameters, rep, times = N * m)
      pes <- pbreplicate(n = N, expr = simulate_pe(parameters = parameters,
                                                   m = m,
                                                   extrapolate_ok = extrapolate_ok,
                                                   shift = shift,
                                                   attempt_fast = attempt_fast), simplify = TRUE) |> as.vector()
      out <- c(spa, list("m" = mid, "x" = nid, "pe" = pes)) |> setDT()
    }
  }
  return(out)
}

#' Simulate Mean Squared Prediction Error for Smoothing Spline
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param N An \code{integer} that represents the number of replicated MSPE calculations
#'          based on each parameter combination.
#' @param m An \code{integer} larger than or equal to \code{1L}.
#'          How many external quality assessment materials should be simulated?
#' @param extrapolate_ok A \code{logical} value. If \code{TRUE}, test data can be
#'                       simulated beyond the empirical range of the clinical sample data.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A}
#'              and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the predictions are attempted speeded up.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'                 If set to \code{TRUE}, simulations are done on multiple cores on
#'                 your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores
#'                  you have on your computer. If you wish to use fewer cores than
#'                  is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here.
#'             If you do not want reproducibility, set this to \code{NULL}.
#' @param percent A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'                If set to \code{TRUE}, all numbers that are possible to present
#'                as percentages are presented as percentages.
#' @param log_mspe A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'                 If set to \code{TRUE}, simulated MSPE values are log-transformed.
#'                 This may be useful for plotting.
#' @param include_ici A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'                    If set to \code{TRUE}, first and third quartiles of simulated
#'                    MSPE values are returned part of the output object.
#' @param trim A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'             If set to \code{TRUE}, trimmed mean with fraction \eqn{0.5\%} is used.
#'
#' @return A \code{data.table}.
#' @export
#'
#' @examples print(1)
simulate_mspe <- function(parameters, N = 100, m = 1, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE, parallel = TRUE, max_cores = 4, seed = 99, percent = TRUE, log_mspe = FALSE, include_ici = FALSE, trim = FALSE){
  pe <- NULL
  multiplier <- ifelse(percent, 100, 1)
  trim_fraction <- ifelse(trim, 0.005, 0)
  simulated_pes <- simulate_pes(parameters = parameters,
                                N = N,
                                m = m,
                                extrapolate_ok = extrapolate_ok,
                                shift = shift,
                                attempt_fast = attempt_fast,
                                parallel = parallel,
                                max_cores = max_cores,
                                seed = seed)
  if(m > 1){
    simulated_pes <- simulated_pes[, list(pe = mean(pe, na.rm = TRUE)),
                                   by = c(names(parameters), "x")]
    if(log_mspe){
      simulated_mspes <- simulated_pes[, list(mspe = mean(log(pe^2), na.rm = TRUE, trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(log(pe^2), probs = 0.25, na.rm = TRUE),
                                              upr = quantile(log(pe^2), probs = 0.75, na.rm = TRUE)),
                                       by = names(parameters)]
    }
    else{
      simulated_mspes <- simulated_pes[, list(mspe = mean(pe^2, na.rm = TRUE, trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(pe^2, probs = 0.25, na.rm = TRUE),
                                              upr = quantile(pe^2, probs = 0.75, na.rm = TRUE)),
                                       by = names(parameters)]
    }

  }
  else{
    if(log_mspe){
      simulated_mspes <- simulated_pes[, list(mspe = mean(log(pe^2), na.rm = TRUE, trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(log(pe^2), probs = 0.25, na.rm = TRUE),
                                              upr = quantile(log(pe^2), probs = 0.75, na.rm = TRUE)),
                                       by = names(parameters)]
    }
    else{
      simulated_mspes <- simulated_pes[, list(mspe = mean(pe^2, na.rm = TRUE, trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(pe^2, probs = 0.25, na.rm = TRUE),
                                              upr = quantile(pe^2, probs = 0.75, na.rm = TRUE)),
                                       by = names(parameters)]
    }
  }

  if(log_mspe){
    #simulated_mspes$mspe <- log(simulated_mspes$mspe)
    #simulated_mspes$lwr <- log(simulated_mspes$lwr)
    #simulated_mspes$upr <- log(simulated_mspes$upr)
  }

  if(!include_ici){
    simulated_mspes$lwr <- simulated_mspes$upr <- NULL
  }

  if(any("cvx" == names(simulated_mspes))){
    simulated_mspes$cvx <- simulated_mspes$cvx * multiplier
  }
  if(any("cvy" == names(simulated_mspes))){
    simulated_mspes$cvy <- simulated_mspes$cvy * multiplier
  }
  if(any("qran" == names(simulated_mspes))){
    simulated_mspes$qran <- simulated_mspes$qran * multiplier
  }

  return(simulated_mspes)

}
