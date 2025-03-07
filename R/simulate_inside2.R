#' @title Simulate Clinical Sample Data and Evaluated Material Data and
#'        Check Whether Inside PI(s)
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame}.
#'                   The parameters and corresponding values to be used in
#'                   the simulation. Either \code{df} or \code{df_max}, must
#'                   be included. See \code{sim_eqa_data()} for other
#'                   parameters.
#' @param m An \code{integer}. Must be larger than or equal to \code{1L}.
#'          The desired number of evaluated materials to simulate. Defaults to
#'          \code{1L}.
#' @param level A \code{double}. Must be between 0 and 1. The desired nominal
#'              confidence level for the approximate prediction intervals.
#'              Defaults to \code{0.95} (\eqn{95\%}).
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
#' @param include_deming A non-missing \code{logical} value. If \code{TRUE},
#'                       the inside prediction interval indicator is also
#'                       calculated for \code{fg} Deming regression. Defaults
#'                       to \code{FALSE}.
#'
#' @details
#' Simulate a clinical sample dataset and an evaluated material dataset. The
#' evaluated material are assumed to be commutable and will thus belong to the
#' same population as the clinical samples. Then the model(s) and corresponding
#' prediction interval(s) is(are) estimated. Finally, one check whether the
#' simulated \code{MP_A} measurement(s) from the evaluated data is inside the
#' prediction interval(s). If it(they) are, \code{inside = 1L}, and otherwise
#' \code{inside = 0L}.
#'
#' If \code{include_deming = TRUE}, these checks are performed for both the
#' Fuller & Gillard Deming regression model and the smoothing spline model.
#'
#' Note: This function is generally not used directly, but as part of the
#' more comprehensive function \code{simulate_confidence_levels2()}.
#'
#' This function will never throw an error. If the invalid arguments are passed
#' the output will instead produce \code{NA}.
#'
#' General note: Any function named \code{simulate_***}, is generally not
#' intended to be used by end-users. They are used in research and validity
#' testing.
#'
#'
#' @return
#' An \code{integer} vector. Is of length \code{2m} if
#' \code{include_deming = TRUE}. Otherwise, length \code{m}.
#'
#' @export
#'
#' @examples
#'
#' # Required packagtes ..
#' library(fasteqa)
#'
#' # Used parameters
#' parameters <- list(n = 30,
#'                    R = 3,
#'                    cvx = 0.01,
#'                    cvy = 0.01,
#'                    df = 6,
#'                    type = 2)
#'
#' # Simulate 100 inside indicator observations (only smoothing splines)
#' ss_pi_insides <- replicate(n = 100,
#'                            expr = simulate_inside2(parameters),
#'                            simplify = TRUE)
#'
#' # Taking the mean produces the approximated confidence level for the
#' # smoothing spline prediction intervals for the given parameters.
#' ss_emp_cl_pi <- mean(ss_pi_insides, na.rm = TRUE)
#'
#' # The result
#' print(ss_emp_cl_pi)
#'

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

  simulated_cs_data <- sim_eqa_data(parameters = parameters,
                                    type = type,
                                    AR = include_deming,
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
    if(include_deming){
      impr <- global_precision_estimates(simulated_cs_data)
      simulated_cs_data <- fun_of_replicates(simulated_cs_data)
      if(eq_parameters$cvx < 1e-8){
        impr$Var_B <- 0
      }
    }
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
    impr <- global_precision_estimates(simulated_cs_data$simulated_data)
    simulated_cs_data$simulated_data <- fun_of_replicates(simulated_cs_data$simulated_data)
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
#' This function checks if evaluated materials are inside the estimated
#' prediction intervals
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
                                                  include_deming = include_deming),
                   simplify = FALSE) |> unlist(use.names = FALSE)
  spa <- lapply(parameters, rep, times = N * m * mul)
  if(include_deming){
    aid <- rep(rep(c("ss", "ev"), each = m), times = N)
    out <- c(spa, list("m" = mid, "x" = nid, "a" = aid, "inside" = ins))
    return(out)
  }
  out <- c(spa, list("m" = mid, "x" = nid, "inside" = ins))
  return(out)
}

#' @title
#' Simulate Clinical Sample Data and Evaluated Material Data and Check Whether
#' Inside PIs (Replicated)
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame}.
#'                   The parameters and corresponding values to be used in
#'                   the simulation. Either \code{df} or \code{df_max}, must
#'                   be included. See \code{sim_eqa_data()} for other
#'                   parameters.
#' @param N An \code{integer}. The desired number of replicated inside checks.
#' @param m An \code{integer}. Must be larger than or equal to \code{1L}.
#'          The desired number of evaluated materials to simulate. Defaults to
#'          \code{1L}.
#' @param level A \code{double}. Must be between 0 and 1. The desired nominal
#'              confidence level for the approximate prediction intervals.
#'              Defaults to \code{0.95} (\eqn{95\%}).
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
#' @param include_deming A non-missing \code{logical} value. If \code{TRUE},
#'                       the inside prediction interval indicator is also
#'                       calculated for \code{fg} Deming regression. Defaults
#'                       to \code{FALSE}.
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
#' Replicates \code{simulate_inside2()} \code{N} times.
#'
#' @details
#' See details in \code{?simulate_inside2()} for the finer details.
#'
#' General note: Any function named \code{simulate_***}, is generally not
#' intended to be used by end-users. They are used in research and validity
#' testing.
#'
#'
#' @return
#' A \code{data.table}. The simulation results for the given \code{parameters}.
#' Expected to have \code{N} \eqn{\times} \code{m} rows. If \code{obs_tau} is
#' given in \code{parameters}, the resulting \code{data.table} is expected to
#' have \code{N} \eqn{\times} \code{m} \eqn{\times} \code{length(obs_tau)} rows.
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
      new_enviroment$sim_eqa_data <- sim_eqa_data
      new_enviroment$global_precision_estimates <- global_precision_estimates
      new_enviroment$fun_of_replicates <- fun_of_replicates
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
                                         "sim_eqa_data",
                                         "global_precision_estimates",
                                         "fun_of_replicates",
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
      new_enviroment$sim_eqa_data <- sim_eqa_data
      new_enviroment$global_precision_estimates <- global_precision_estimates
      new_enviroment$fun_of_replicates <- fun_of_replicates
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
                                         "sim_eqa_data",
                                         "global_precision_estimates",
                                         "fun_of_replicates",
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







