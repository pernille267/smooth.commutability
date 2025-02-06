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
#' @param include_deming A \code{logical} value. If \code{TRUE}, the inside rates based on Deming regression are also included.
#'
#' @return An \code{integer} vector of length \code{m}.
#' @export
#'
#' @examples print(1)
simulate_inside2 <- function(parameters,
                             m = 1,
                             level = 0.95,
                             extrapolate_ok = FALSE,
                             shift = FALSE,
                             attempt_fast = FALSE,
                             include_deming = FALSE){
  impr <- NULL
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

  simulated_cs_data <- simulate_eqa_data2(parameters = parameters,
                                          type = type,
                                          AR = include_deming,
                                          include_parameters = TRUE,
                                          shift = shift_roles)
  if(!is.null(obs_tau)){
    eq_parameters <- simulated_cs_data$parameters
    eq_parameters$n <- m
    eq_parameters$obs_tau <- obs_tau
    simulated_cs_data <- simulated_cs_data$simulated_data
    if(include_deming){
      impr <- global_precision_estimates2(simulated_cs_data)
      simulated_cs_data <- fun_of_replicates2(simulated_cs_data)
      if(eq_parameters$cvx < 1e-8){
        impr$Var_B <- 0
      }

    }
    simulated_eq_data <- simulate_eqa_data2(parameters = eq_parameters,
                                            type = type,
                                            AR = FALSE,
                                            include_parameters = FALSE,
                                            shift = shift_roles)
    output_ss <- tryCatch(expr = predict_smoothing_spline(data = simulated_cs_data,
                                                          new_data = simulated_eq_data,
                                                          weighted = FALSE,
                                                          df = if(any("df" == names(parameters))){parameters$df}else{NULL},
                                                          lambda = parameters$lambda,
                                                          df_max = parameters$df_max,
                                                          R_ratio = 1,
                                                          level = level,
                                                          negative_ok = TRUE,
                                                          attempt_fast = attempt_fast,
                                                          rounding = 3L,
                                                          additional_parameters = additional_parameters)$inside,
                          error = function(e) rep(NA_integer_, m))
    if(include_deming){
      output_ev <- inside_deming2(data = simulated_cs_data,
                                  new_data = simulated_eq_data,
                                  imprecision_estimates = impr,
                                  R = parameters$R,
                                  R_ratio = 1,
                                  level = level)
      return(c(output_ss, output_ev))
    }
    return(output_ss)
  }
  eq_parameters <- simulated_cs_data$parameters
  eq_parameters$n <- m
  if(include_deming){
    impr <- global_precision_estimates2(simulated_cs_data$simulated_data)
    simulated_cs_data$simulated_data <- fun_of_replicates2(simulated_cs_data$simulated_data)
    if(eq_parameters$cvx < 1e-8){
      impr$Var_B <- 0
    }
  }
  if(!extrapolate_ok){
    eq_parameters$cil <- min(simulated_cs_data$simulated_data$MP_B, na.rm = TRUE)
    eq_parameters$ciu <- max(simulated_cs_data$simulated_data$MP_B, na.rm = TRUE)
    if(eq_parameters$ciu <= eq_parameters$cil){
      if(include_deming){
        return(rep(NA_integer_, m * 2))
      }
      else{
        return(rep(NA_integer_, m))
      }

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
                                                        level = level,
                                                        negative_ok = TRUE,
                                                        attempt_fast = attempt_fast,
                                                        rounding = 3L,
                                                        additional_parameters = additional_parameters)$inside,
                        error = function(e) rep(NA_integer_, m))
  if(include_deming){
    output_ev <- inside_deming2(data = simulated_cs_data,
                                new_data = simulated_eq_data,
                                imprecision_estimates = impr,
                                R = parameters$R,
                                R_ratio = 1,
                                level = level)
    return(c(output_ss, output_ev))
  }
  return(output_ss)
}

#' Internal Function for Inside Calculation
#'
#' This function checks if external quality assessment materials are inside estimated prediction intervals
#' @keywords internal
simulate_insides2_internal <- function(parameters,
                                       N = 100,
                                       m = 1,
                                       level = 0.95,
                                       extrapolate_ok = FALSE,
                                       shift = FALSE,
                                       attempt_fast = FALSE,
                                       include_deming = FALSE){
  mul <- if(include_deming){2}else{1}
  mid <- rep(1:m, times = N * mul)
  nid <- rep(1:N, each = m * mul)
  ins <- replicate(n = N, expr = simulate_inside2(parameters = parameters,
                                                  m = m,
                                                  level = level,
                                                  extrapolate_ok = extrapolate_ok,
                                                  shift = shift,
                                                  attempt_fast = attempt_fast,
                                                  include_deming = include_deming), simplify = FALSE) |> unlist(use.names = FALSE)
  spa <- lapply(parameters, rep, times = N * m * mul)
  if(include_deming){
    aid <- rep(rep(c("ss", "ev"), each = m), times = N)
    out <- c(spa, list("m" = mid, "x" = nid, "a" = aid, "inside" = ins))
    return(out)
  }
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
#' @param include_deming A \code{logical} value. If \code{TRUE}, the inside rates based on Deming regression are also included.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here. If you do not want reproducibility, set this to \code{NULL}.
#'
#' @return A \code{list}.
#' @export
#'
#' @examples print(1)
simulate_insides2 <- function(parameters,
                              N = 100,
                              m = 1,
                              level = 0.95,
                              extrapolate_ok = FALSE,
                              shift = FALSE,
                              attempt_fast = FALSE,
                              include_deming = FALSE,
                              parallel = TRUE,
                              max_cores = 4,
                              seed = 99){
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
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$include_deming <- include_deming
      new_enviroment$simulate_insides2_internal <- simulate_insides2_internal
      new_enviroment$simulate_inside2 <- simulate_inside2
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2
      new_enviroment$global_precision_estimates2 <- global_precision_estimates2
      new_enviroment$fun_of_replicates2 <- fun_of_replicates2
      new_enviroment$inside_deming2 <- inside_deming2


      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "m",
                                         "parameters_list",
                                         "level",
                                         "extrapolate_ok",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "include_deming",
                                         "simulate_eqa_data2",
                                         "global_precision_estimates2",
                                         "fun_of_replicates2",
                                         "inside_deming2",
                                         "simulate_insides2_internal",
                                         "simulate_inside2"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_insides2_internal(parameters = x,
                                                                                        N = N,
                                                                                        m = m,
                                                                                        level = level,
                                                                                        extrapolate_ok = extrapolate_ok,
                                                                                        shift = shift,
                                                                                        attempt_fast = attempt_fast,
                                                                                        include_deming = include_deming), cl = cl)
      out <- rbindlist(out)
    }
    else if(!isTRUE(parallel)){
      out <- pblapply(X = parameters_list, FUN = function(x) simulate_insides2_internal(parameters = x,
                                                                                        N = N,
                                                                                        m = m,
                                                                                        level = level,
                                                                                        extrapolate_ok = extrapolate_ok,
                                                                                        shift = shift,
                                                                                        attempt_fast = attempt_fast,
                                                                                        include_deming = include_deming))
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
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$include_deming <- include_deming
      new_enviroment$simulate_insides2_internal <- simulate_insides2_internal
      new_enviroment$simulate_inside2 <- simulate_inside2
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2
      new_enviroment$global_precision_estimates2 <- global_precision_estimates2
      new_enviroment$fun_of_replicates2 <- fun_of_replicates2
      new_enviroment$inside_deming2 <- inside_deming2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "m",
                                         "parameters",
                                         "level",
                                         "extrapolate_ok",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "include_deming",
                                         "simulate_insides2_internal",
                                         "simulate_inside2",
                                         "simulate_eqa_data2",
                                         "global_precision_estimates2",
                                         "fun_of_replicates2",
                                         "inside_deming2"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      mul <- if(include_deming){2}else{1}
      mid <- rep(1:m, times = N * mul)
      nid <- rep(1:N, each = m * mul)
      spa <- lapply(parameters, rep, times = N * m * mul)
      ins <- pbreplicate(n = N, expr = simulate_inside2(parameters = parameters,
                                                        m = m,
                                                        level = level,
                                                        extrapolate_ok = extrapolate_ok,
                                                        shift = shift,
                                                        attempt_fast = attempt_fast,
                                                        include_deming = include_deming), cl = cl, simplify = FALSE) |> unlist(use.names = FALSE)
      if(include_deming){
        aid <- rep(rep(c("ss", "ev"), each = m), times = N)
        out <- c(spa, list("m" = mid, "x" = nid, "a" = aid, "inside" = ins)) |> setDT()
        return(out)
      }
      out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins)) |> setDT()

    }
    else if(!isTRUE(parallel)){
      mul <- if(include_deming){2}else{1}
      mid <- rep(1:m, times = N * mul)
      nid <- rep(1:N, each = m * mul)
      spa <- lapply(parameters, rep, times = N * m * mul)
      ins <- pbreplicate(n = N, expr = simulate_inside2(parameters = parameters,
                                                        m = m,
                                                        level = level,
                                                        extrapolate_ok = extrapolate_ok,
                                                        shift = shift,
                                                        attempt_fast = attempt_fast,
                                                        include_deming = include_deming), simplify = FALSE) |> unlist(use.names = FALSE)
      if(include_deming){
        aid <- rep(rep(c("ss", "ev"), each = m), times = N)
        out <- c(spa, list("m" = mid, "x" = nid, "a" = aid, "inside" = ins)) |> setDT()
        return(out)
      }
      out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins)) |> setDT()
    }
  }
  return(out)
}







