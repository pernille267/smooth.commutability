#' @title
#' Simulate Smoothing Spline Prediction Errors
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame}.
#'                   The parameters and corresponding values to be used in
#'                   the simulation. Either \code{df} or \code{df_max}, must
#'                   be included. See \code{sim_eqa_data()} for other
#'                   parameters.
#' @param m An \code{integer}. Must be larger than or equal to \code{1L}.
#' @param extrapolate_ok A non-missing \code{logical} value. If \code{TRUE},
#'                       simulated evaluated material data are allowed to
#'                       produce measurements beyond the support of the
#'                       clinical samples IVD-MD measurements. Defaults to
#'                       \code{FALSE}.
#' @param shift A non-missing \code{logical} value. If \code{TRUE}, the roles
#'              of \code{MP_A} and \code{MP_B} are allowed to shift. Can only
#'              happen if \code{lambda < 0.5}. Defaults to \code{FALSE}.
#' @param attempt_fast A non-missing \code{logical} value. If \code{TRUE}, the
#'                     \code{smooth.spline()} is used to make predictions.
#'
#' @description
#' Simulates smoothing spline prediction errors for each of the \code{m}
#' evaluated materials.
#'
#' @details
#' Simulate a clinical sample dataset and an evaluated material dataset. The
#' evaluated material are assumed to be commutable and will thus belong to the
#' same population as the clinical samples. Then the corresponding
#' prediction error(s) is(are) calculated. The prediction errors is of course
#' calculated by the arithmetic difference between the observed value(s) of
#' the evaluated material and the predicted value using the smoothing spline
#' model.
#'
#' Note: This function is generally not used directly, but as part of the
#' more comprehensive function \code{simulate_mspe()}.
#'
#' This function will never throw an error. If the invalid arguments are passed
#' the output will instead produce \code{NA}.
#'
#' General note: Any function named \code{simulate_***}, is generally not
#' intended to be used by end-users. They are used in research and validity
#' testing.
#'
#' @return
#' A \code{numeric} vector of length \code{m}. The simulated prediction errors.
#' @export
#'
#' @examples
#' # Required packages
#' library(fasteqa)
#'
#' # Used parameters
#' parameters <- list(n = 25,
#'                    R = 3,
#'                    cvx = 0.01,
#'                    cvy = 0.01,
#'                    obs_tau = 6,
#'                    df = 6,
#'                    type = 2)
#'
#' # Simulate 1000 smoothing spline prediction errors at y = 6
#' ss_pes <- replicate(n = 1000,
#'                     expr = simulate_pe(parameters),
#'                     simplify = TRUE)
#'
#' # Calculate prediction bias at y = 6
#' print(round(mean(ss_pes, na.rm = TRUE), 4L))
#'
#' # Calculate variance of prediction error at y = 6
#' print(round(var(ss_pes, na.rm = TRUE), 4L))
#'
#'
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
  simulated_cs_data <- sim_eqa_data(parameters = parameters,
                                    type = type,
                                    AR = FALSE,
                                    include_parameters = TRUE)

  if(shift_roles){
    temp_MP_A <- simulated_cs_data$simulated_data$MP_A
    simulated_cs_data$simulated_data$MP_A <- simulated_cs_data$simulated_data$MP_B
    simulated_cs_data$simulated_data$MP_B <- temp_MP_A
  }


  if(!is.null(obs_tau)){
    eq_parameters <- simulated_cs_data$parameters
    eq_parameters$n <- m
    eq_parameters$obs_tau <- obs_tau
    simulated_cs_data <- simulated_cs_data$simulated_data
    simulated_eq_data <- sim_eqa_data(parameters = eq_parameters,
                                      type = type,
                                      AR = FALSE,
                                      include_parameters = FALSE)

    if(shift_roles){
      temp_MP_A <- simulated_eq_data$MP_A
      simulated_eq_data$MP_A <- simulated_eq_data$MP_B
      simulated_eq_data$MP_B <- temp_MP_A
    }

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
  simulated_eq_data <- sim_eqa_data(parameters = eq_parameters,
                                    type = type,
                                    AR = FALSE,
                                    include_parameters = FALSE)

  if(shift_roles){
    temp_MP_A <- simulated_eq_data$MP_A
    simulated_eq_data$MP_A <- simulated_eq_data$MP_B
    simulated_eq_data$MP_B <- temp_MP_A
  }

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
#' This function replicates prediction error calculation \code{N} times for
#' \code{parameters}
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

#' @title
#' Replicate Calculation of Prediction Error for Each Parameter Combination
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame}.
#'                   The parameters and corresponding values to be used in
#'                   the simulation. Either \code{df} or \code{df_max}, must
#'                   be included. See \code{sim_eqa_data()} for other
#'                   parameters.
#' @param N An \code{integer}. The desired number of replicated prediction error
#'          calculations for each parameter combination in \code{parameters}.
#' @param m An \code{integer}. Must be larger than or equal to \code{1L}.
#'          The desired number of evaluated materials to simulate. Defaults to
#'          \code{1L}.
#' @param extrapolate_ok A non-missing \code{logical} value. If \code{TRUE},
#'                       simulated evaluated material data are allowed to
#'                       produce measurements beyond the support of the
#'                       clinical samples IVD-MD measurements. Defaults to
#'                       \code{FALSE}.
#' @param shift A non-missing \code{logical} value. If \code{TRUE}, the roles
#'              of \code{MP_A} and \code{MP_B} are allowed to shift. Can only
#'              happen if \code{lambda < 0.5}. Defaults to \code{FALSE}.
#' @param attempt_fast A non-missing \code{logical} value. If \code{TRUE}, the
#'                     \code{smooth.spline()} is used to make predictions.
#' @param parallel A non-missing \code{logical} value. If \code{TRUE},
#'                 simulations are done in parallel. May reduce runtime.
#' @param max_cores An \code{integer}. Must be between \code{2} and the number
#'                  of cores you have on your computer. If fewer cores than
#'                  the total number of cores on your computer is desired,
#'                  specify this here.
#' @param seed An \code{integer}. Ensuring reproducibility even if
#'             \code{parallel = TRUE}. If reproducibility is unimportant, set
#'             this to \code{NULL}.
#'
#' @description
#' Replicates \code{simulate_pe()} \code{N} times.
#'
#' @details
#' Additional details...
#'
#'
#' @return
#' A \code{data.table}. The approximated mean squared prediction errors (MSPEs)
#' for each unique combination of parameters in \code{parameters}. The
#' resulting \code{data.table} is expected to have
#' \code{N} \eqn{\times} \code{m} \eqn{\times} \code{length(obs_tau)} rows.
#' @export
#'
#' @examples
#' # Required packages
#' library(fasteqa)
#'
#' # Used parameters
#' parameters <- list(n = 25,
#'                    R = 3,
#'                    cvx = 0.01,
#'                    cvy = 0.01,
#'                    cil = 2,
#'                    ciu = 10,
#'                    obs_tau = 9,
#'                    type = 3,
#'                    df = 5)
#'
#' # Simulated smoothing spline prediction errors
#' ss_pes <- simulate_pes(parameters = parameters,
#'                        N = 1000,
#'                        parallel = FALSE)
#'
#' # Approximate prediction bias at y = 9
#' print(round(mean(ss_pes$pe, na.rm = TRUE), 4L))
#'
#' # Approximate MSPE at y = 9
#' print(round(mean(ss_pes$pe^2, na.rm = TRUE), 4L))
#'

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
      new_enviroment$sim_eqa_data <- sim_eqa_data

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "m",
                                         "parameters_list",
                                         "extrapolate_ok",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "sim_eqa_data",
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
      new_enviroment$sim_eqa_data <- sim_eqa_data

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
                                         "sim_eqa_data"), envir = new_enviroment)
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

#' @title
#' Simulate Smoothing Spline Mean Squared Prediction Errors
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame}.
#'                   The parameters and corresponding values to be used in
#'                   the simulation. Either \code{df} or \code{df_max}, must
#'                   be included. See \code{sim_eqa_data()} for other
#'                   parameters.
#' @param N An \code{integer}. The desired number of replicated prediction error
#'          calculations for each parameter combination in \code{parameters}.
#' @param m An \code{integer}. Must be larger than or equal to \code{1L}.
#'          The desired number of evaluated materials to simulate. Defaults to
#'          \code{1L}.
#' @param extrapolate_ok A non-missing \code{logical} value. If \code{TRUE},
#'                       simulated evaluated material data are allowed to
#'                       produce measurements beyond the support of the
#'                       clinical samples IVD-MD measurements. Defaults to
#'                       \code{FALSE}.
#' @param shift A non-missing \code{logical} value. If \code{TRUE}, the roles
#'              of \code{MP_A} and \code{MP_B} are allowed to shift. Can only
#'              happen if \code{lambda < 0.5}. Defaults to \code{FALSE}.
#' @param attempt_fast A non-missing \code{logical} value. If \code{TRUE}, the
#'                     \code{smooth.spline()} is used to make predictions.
#'                     Defaults to \code{FALSE}.
#' @param parallel A non-missing \code{logical} value. If \code{TRUE},
#'                 simulations are done in parallel. May reduce runtime.
#'                 Defaults to \code{FALSE}.
#' @param max_cores An \code{integer}. Must be between \code{2} and the number
#'                  of cores you have on your computer. If fewer cores than
#'                  the total number of cores on your computer is desired,
#'                  specify this here. Defaults to \code{4L}.
#' @param seed An \code{integer}. Ensuring reproducibility even if
#'             \code{parallel = TRUE}. If reproducibility is unimportant, set
#'             this to \code{NULL}. Defaults to \code{99L}.
#' @param percent A non-missing \code{logical} value. If \code{TRUE}, all
#'                results possible to present as percentages are presented as
#'                percentages. Defaults to \code{TRUE}.
#' @param log_mspe A non-missing \code{logical} value. If \code{TRUE},
#'                 simulated MSPE values are log-transformed. May be useful
#'                 for plotting purposes if raw MSPE values near \code{0}.
#'                 Defaults to \code{FALSE}.
#' @param include_ici A non-missing \code{logical} value. If \code{TRUE}, the
#'                    first and third quartile of simulated MSPE values
#'                    included in the output. Defaults to \code{FALSE}.
#' @param trim A non-missing \code{logical}. If \code{TRUE}, MSPE values are
#'             calculated after trimming extreme and potentially errenous
#'             prediction errors. Defaults to \code{FALSE}.
#'
#' @description
#' Simulates MSPE values for each parameter combination in \code{parameters}.
#'
#' @details
#' Replicates \code{simulate_pe()} \code{N} times for each parameter
#' combination in \code{parameters}. For each combination, the arithmetic mean
#' is applied to the squared simulated prediction errors, and the end result is
#' the simulated MSPE values.
#'
#' General note: Any function named \code{simulate_***}, is generally not
#' intended to be used by end-users. They are used in research and validity
#' testing.
#'
#'
#' @return
#' A \code{data.table}. The resulting simulated MSPE values for every unique
#' parameter combination in \code{parameters}. If \code{include_ici = TRUE},
#' quartiles of the simulated squared prediction errors are included in the
#' output and named \code{lwr} and \code{upr}, respectively.
#' @export
#'
#' @examples
#' # Required packages
#' library(fasteqa)
#' library(data.table)
#'
#' # Simulation parameters
#' n <- c(25, 40)
#' R <- c(3, 4)
#' cvx <- 0.01
#' cvy <- 0.01
#' type <- c(1, 2, 3)
#' df_max <- 7.5
#' tau_obs <- c(3, 6, 9)
#' parameters <- CJ(n, R, cvx, cvy, df_max, type, tau_obs)
#'
#' # Number of unique parameter combinations
#' print(nrow(parameters))
#'
#' # Approximate smoothing spline MSPE values (using trimming)
#' # Takes approximately 1.5 minutes to run
#' output <- simulate_mspe(parameters = parameters,
#'                         N = 1e3,
#'                         m = 1,
#'                         parallel = FALSE,
#'                         percent = TRUE,
#'                         trim = TRUE)
#'
#' # Rounding results for cleaner output
#' output$mspe <- round(output$mspe, 4L)
#'
#' # The output
#' print(output)
#'
simulate_mspe <- function(parameters,
                          N = 100,
                          m = 1,
                          extrapolate_ok = FALSE,
                          shift = FALSE,
                          attempt_fast = FALSE,
                          parallel = FALSE,
                          max_cores = 4,
                          seed = 99,
                          percent = TRUE,
                          log_mspe = FALSE,
                          include_ici = FALSE,
                          trim = FALSE){
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
      simulated_mspes <- simulated_pes[, list(mspe = mean(log(pe^2),
                                                          na.rm = TRUE,
                                                          trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(log(pe^2),
                                                             probs = 0.25,
                                                             na.rm = TRUE),
                                              upr = quantile(log(pe^2),
                                                             probs = 0.75,
                                                             na.rm = TRUE)),
                                       by = names(parameters)]
    }
    else{
      simulated_mspes <- simulated_pes[, list(mspe = mean(pe^2,
                                                          na.rm = TRUE,
                                                          trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(pe^2,
                                                             probs = 0.25,
                                                             na.rm = TRUE),
                                              upr = quantile(pe^2,
                                                             probs = 0.75,
                                                             na.rm = TRUE)),
                                       by = names(parameters)]
    }

  }
  else{
    if(log_mspe){
      simulated_mspes <- simulated_pes[, list(mspe = mean(log(pe^2),
                                                          na.rm = TRUE,
                                                          trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(log(pe^2),
                                                             probs = 0.25,
                                                             na.rm = TRUE),
                                              upr = quantile(log(pe^2),
                                                             probs = 0.75,
                                                             na.rm = TRUE)),
                                       by = names(parameters)]
    }
    else{
      simulated_mspes <- simulated_pes[, list(mspe = mean(pe^2,
                                                          na.rm = TRUE,
                                                          trim = trim_fraction),
                                              number_of_observations = sum(!is.na(pe)),
                                              lwr = quantile(pe^2,
                                                             probs = 0.25,
                                                             na.rm = TRUE),
                                              upr = quantile(pe^2,
                                                             probs = 0.75,
                                                             na.rm = TRUE)),
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
