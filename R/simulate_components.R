#' Simulate Smoothing Spline Components
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the smoothing spline fitting and predictions are attempted speeded up.
#' @param include_pe A \code{logical} value. If \code{TRUE}, predictions and prediction errors are returned.
#'
#' @return A named \code{integer} vector containing the relevant simulated components based on \code{parameters}
#' @export
#'
#' @examples print(1)
simulate_components <- function(parameters,
                                shift = FALSE,
                                attempt_fast = FALSE,
                                include_pe = FALSE){
  m <- 1
  type <- NULL
  obs_tau <- NULL
  additional_parameters <- list("method" = "loocv",
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
  output_ss <- tryCatch(expr = smoothing_spline(data = simulated_cs_data$simulated_data,
                                                weights = 1,
                                                df = if(any("df" == names(parameters))){parameters$df}else{NULL},
                                                lambda = parameters$lambda,
                                                df_max = parameters$df_max,
                                                attempt_fast = attempt_fast,
                                                na_rm = TRUE,
                                                additional_parameters = additional_parameters),
                        error = function(e) NULL)
  output_var_eps <- NA_real_
  output_df <- NA_real_
  if(!is.null(output_ss)){
    output_var_eps <- output_ss$var_eps
    output_df <- output_ss$df
    if(output_df < 2 - 1e-6 | output_df > length(output_ss$x)){
      output_var_eps <- NA_real_
      output_df <- NA_real_
    }
  }


  if(!is.null(obs_tau)){
    if(is.character(obs_tau)){
      obs_tau <- runif(1, simulated_cs_data$parameters$cil, simulated_cs_data$parameters$ciu)
    }
    m <- length(obs_tau)
    eq_parameters <- simulated_cs_data$parameters
    eq_parameters$n <- m
    eq_parameters$obs_tau <- obs_tau
    simulated_cs_data <- simulated_cs_data$simulated_data

    simulated_eq_data <- simulate_eqa_data2(parameters = eq_parameters, type = type, AR = FALSE, include_parameters = FALSE, shift = shift_roles)
    output_pr <- tryCatch(expr = predict_smoothing_spline(data = simulated_cs_data,
                                                          new_data = simulated_eq_data,
                                                          weighted = FALSE,
                                                          df = if(any("df" == names(parameters))){parameters$df}else{NULL},
                                                          lambda = parameters$lambda,
                                                          df_max = parameters$df_max,
                                                          R_ratio = 1,
                                                          negative_ok = TRUE,
                                                          attempt_fast = attempt_fast,
                                                          include_prediction_variance = TRUE,
                                                          rounding = 3L,
                                                          additional_parameters = additional_parameters),
                          error = function(e) NULL)
    vpe <- rep(NA_real_, m)
    pe <- rep(NA_real_, m)
    pr <- rep(NA_real_, m)
    if(!is.null(output_pr)){
      vpe <- output_pr$var_pred_error
      pr <- output_pr$prediction
      pe <- simulated_eq_data$MP_A - output_pr$prediction
    }
    var_pred <- vpe - output_var_eps
    if(m == 1){
      if(include_pe){
        return(c("pred_error" = pe, "pred" = pr, "var_eps" = output_var_eps, "var_pred" = var_pred, "var_pred_error" = vpe, "df" = output_df, "t" = pe / sqrt(vpe)))
      }
      else{
        return(c("var_eps" = output_var_eps, "var_pred" = var_pred, "var_pred_error" = vpe, "df" = output_df, "t" = pe / sqrt(vpe)))
      }

    }
    else{
      if(include_pe){
        return(c(pe, pr, var_pred, vpe, pe / sqrt(vpe)))
      }
      else{
        return(c(var_pred, vpe, pe / sqrt(vpe)))
      }

    }
  }
  return(c("var_eps" = output_var_eps, "df" = output_df))
}

#' Internal Function for Component Simulation
#'
#' This function repeats the simulate_components() function \code{N} times for a single set of simulation parameters.
#' @keywords internal
replicate_simulate_components_internal <- function(parameters,
                                                   N = 100,
                                                   shift = FALSE,
                                                   attempt_fast = FALSE,
                                                   include_pe = FALSE){
  comps_test <- simulate_components(parameters = parameters,
                                    shift = shift,
                                    attempt_fast = attempt_fast,
                                    include_pe = include_pe)
  if(!is.null(names(comps_test))){
    mid <- rep(1, N * length(comps_test))
    nid <- rep(1:N, each = length(comps_test))
    aid <- rep(names(comps_test), N)
    comps <- replicate(n = N, expr = unname(simulate_components(parameters = parameters,
                                                                shift = shift,
                                                                attempt_fast = attempt_fast,
                                                                include_pe = include_pe)), simplify = FALSE) |> unlist(use.names = FALSE)
    spa <- lapply(parameters, rep, times = N * length(comps_test))
    out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps))
    return(out)
  }
  if(include_pe){
    mid <- rep(1:(length(comps_test) / 5), N * 5)
    nid <- rep(1:N, each = length(comps_test))
    aid <- rep(rep(c("pred_error", "pred", "var_pred", "var_pred_error", "t"), each = length(comps_test) / 5), N)
    comps <- replicate(n = N, expr = unname(simulate_components(parameters = parameters,
                                                                shift = shift,
                                                                attempt_fast = attempt_fast,
                                                                include_pe = include_pe)), simplify = FALSE) |> unlist(use.names = FALSE)
    spa <- lapply(parameters, rep, times = N * length(comps_test))
    out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps))
    return(out)
  }
  mid <- rep(1:(length(comps_test) / 3), N * 3)
  nid <- rep(1:N, each = length(comps_test))
  aid <- rep(rep(c("var_pred", "var_pred_error", "t"), each = length(comps_test) / 3), N)
  comps <- replicate(n = N, expr = unname(simulate_components(parameters = parameters,
                                                              shift = shift,
                                                              attempt_fast = attempt_fast,
                                                              include_pe = include_pe)), simplify = FALSE) |> unlist(use.names = FALSE)
  spa <- lapply(parameters, rep, times = N * length(comps_test))
  out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps))
  return(out)

}

#' Replicate Simulation of Smoothing Spline Components
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param N An \code{integer} that represents the number of replicated component estimates.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the calculations are attempted speeded up by using some of the \code{smooth.spline()} functionality.
#' @param include_pe A \code{logical} value. If \code{TRUE}, predictions and prediction errors are returned.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here. If you do not want reproducibility, set this to \code{NULL}.
#'
#' @return A \code{data.table} object, grouped by the unique parameter combinations as well as component types
#' @export
#'
#' @examples print(1)
replicate_simulate_components <- function(parameters,
                                          N = 100,
                                          shift = FALSE,
                                          attempt_fast = FALSE,
                                          include_pe = FALSE,
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
      new_enviroment$parameters_list <- parameters_list
      new_enviroment$shift <- shift
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$include_pe <- include_pe
      new_enviroment$replicate_simulate_components_internal <- replicate_simulate_components_internal
      new_enviroment$simulate_components <- simulate_components
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "parameters_list",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "include_pe",
                                         "simulate_eqa_data2",
                                         "simulate_components",
                                         "replicate_simulate_components_internal"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      out <- pblapply(X = parameters_list, FUN = function(x) replicate_simulate_components_internal(parameters = x,
                                                                                                    N = N,
                                                                                                    shift = shift,
                                                                                                    attempt_fast = attempt_fast,
                                                                                                    include_pe = include_pe), cl = cl)
      out <- rbindlist(out)
    }
    else if(!isTRUE(parallel)){
      out <- pblapply(X = parameters_list, FUN = function(x) replicate_simulate_components_internal(parameters = x,
                                                                                                    N = N,
                                                                                                    shift = shift,
                                                                                                    attempt_fast = attempt_fast,
                                                                                                    include_pe = include_pe))
      out <- rbindlist(out)
    }
  }
  else if(!multiple_pars){
    if(isTRUE(parallel)){
      new_enviroment <- new.env()
      new_enviroment$N <- N
      new_enviroment$parameters <- parameters
      new_enviroment$shift <- shift
      new_enviroment$seed <- seed
      new_enviroment$attempt_fast <- attempt_fast
      new_enviroment$include_pe <- include_pe
      new_enviroment$simulate_components <- simulate_components
      new_enviroment$simulate_eqa_data2 <- simulate_eqa_data2

      cl <- makeCluster(min(detectCores() - 1, max_cores))
      on.exit(stopCluster(cl))
      clusterExport(cl = cl, varlist = c("N",
                                         "parameters",
                                         "shift",
                                         "seed",
                                         "attempt_fast",
                                         "include_pe",
                                         "simulate_eqa_data2",
                                         "simulate_components"), envir = new_enviroment)
      if(is.null(seed)){
        clusterEvalQ(cl = cl, expr = {library(data.table)})
      }
      else{
        clusterEvalQ(cl = cl, expr = {library(data.table);set.seed(seed)})
      }
      comps_test <- simulate_components(parameters = parameters,
                                        shift = shift,
                                        attempt_fast = attempt_fast,
                                        include_pe = include_pe)
      if(!is.null(names(comps_test))){
        mid <- rep(1, N * length(comps_test))
        nid <- rep(1:N, each = length(comps_test))
        aid <- rep(names(comps_test), N)
        comps <- pbreplicate(n = N, expr = unname(simulate_components(parameters = parameters,
                                                                      shift = shift,
                                                                      attempt_fast = attempt_fast,
                                                                      include_pe = include_pe)), simplify = FALSE, cl = cl) |> unlist(use.names = FALSE)
        spa <- lapply(parameters, rep, times = N * length(comps_test))
        out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
        return(out)
      }
      if(include_pe){
        mid <- rep(1:(length(comps_test) / 5), N * 5)
        nid <- rep(1:N, each = length(comps_test))
        aid <- rep(rep(c("pred_error", "pred", "var_pred", "var_pred_error", "t"), each = length(comps_test) / 5), N)
        comps <- pbreplicate(n = N, expr = simulate_components(parameters = parameters,
                                                               shift = shift,
                                                               attempt_fast = attempt_fast,
                                                               include_pe = include_pe), simplify = FALSE, cl = cl) |> unlist(use.names = FALSE)
        spa <- lapply(parameters, rep, times = N * length(comps_test))
        if(!is.null(spa$obs_tau)){
          spa$obs_tau <- spa$obs_tau[1:(length(spa[[1]]))]
        }
        out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
        return(out)
      }
      mid <- rep(1:(length(comps_test) / 3), N * 3)
      nid <- rep(1:N, each = length(comps_test))
      aid <- rep(rep(c("var_pred", "var_pred_error", "t"), each = length(comps_test) / 3), N)
      comps <- pbreplicate(n = N, expr = simulate_components(parameters = parameters,
                                                             shift = shift,
                                                             attempt_fast = attempt_fast,
                                                             include_pe = include_pe), simplify = FALSE, cl = cl) |> unlist(use.names = FALSE)
      spa <- lapply(parameters, rep, times = N * length(comps_test))
      if(!is.null(spa$obs_tau)){
        spa$obs_tau <- spa$obs_tau[1:(length(spa[[1]]))]
      }
      out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
      return(out)

    }
    else if(!isTRUE(parallel)){
      comps_test <- simulate_components(parameters = parameters,
                                        shift = shift,
                                        attempt_fast = attempt_fast,
                                        include_pe = include_pe)
      if(!is.null(names(comps_test))){
        mid <- rep(1, N * length(comps_test))
        nid <- rep(1:N, each = length(comps_test))
        aid <- rep(names(comps_test), N)
        comps <- pbreplicate(n = N, expr = unname(simulate_components(parameters = parameters,
                                                                      shift = shift,
                                                                      attempt_fast = attempt_fast,
                                                                      include_pe = include_pe)), simplify = FALSE) |> unlist(use.names = FALSE)
        spa <- lapply(parameters, rep, times = N * length(comps_test))
        out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
        return(out)
      }
      if(include_pe){
        mid <- rep(1:(length(comps_test) / 5), N * 5)
        nid <- rep(1:N, each = length(comps_test))
        aid <- rep(rep(c("pred_error", "pred", "var_pred", "var_pred_error", "t"), each = length(comps_test) / 5), N)
        comps <- pbreplicate(n = N, expr = simulate_components(parameters = parameters,
                                                               shift = shift,
                                                               attempt_fast = attempt_fast,
                                                               include_pe = include_pe), simplify = FALSE) |> unlist(use.names = FALSE)
        spa <- lapply(parameters, rep, times = N * length(comps_test))
        if(!is.null(spa$obs_tau)){
          spa$obs_tau <- spa$obs_tau[1:(length(spa[[1]]))]
        }
        out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
        return(out)

      }
      mid <- rep(1:(length(comps_test) / 3), N * 3)
      nid <- rep(1:N, each = length(comps_test))
      aid <- rep(rep(c("var_pred", "var_pred_error", "t"), each = length(comps_test) / 3), N)
      comps <- pbreplicate(n = N, expr = simulate_components(parameters = parameters,
                                                             shift = shift,
                                                             attempt_fast = attempt_fast,
                                                             include_pe = include_pe), simplify = FALSE) |> unlist(use.names = FALSE)
      spa <- lapply(parameters, rep, times = N * length(comps_test))
      if(!is.null(spa$obs_tau)){
        spa$obs_tau <- spa$obs_tau[1:(length(spa[[1]]))]
      }
      out <- c(spa, list("m" = mid, "x" = nid, "component" = aid, "value" = comps)) |> setDT()
      return(out)
    }
  }
  return(out)
}

#' Simulation of Statistics of The Smoothing Spline Components
#'
#' @param parameters A \code{list}, \code{data.table} or \code{data.frame} with column names representing the parameters to be included in the simulations.
#'                   Must in general include \code{n}, \code{R}, \code{cvx}, \code{cvy}, and \code{df} or \code{df_max}. If \code{type} is included, data will
#'                   be simulated from custom non-linear functions. Note that \code{cil} and \code{ciu} also must be included if \code{type} is included.
#' @param N An \code{integer} that represents the number of replicated component estimates.
#' @param shift A \code{logical} value. If \code{TRUE}, the roles of \code{MP_A} and \code{MP_B} shift if \code{lambda < 0.5}.
#' @param attempt_fast A \code{logical} value. If \code{TRUE}, the calculations are attempted speeded up by using some of the \code{smooth.spline()} functionality.
#' @param include_pe A \code{logical} value. If \code{TRUE}, predictions and prediction errors are returned.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here. If you do not want reproducibility, set this to \code{NULL}.
#' @param percent A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, all numbers that are possible to present as percentages are presented as percentages.
#' @param remove_extremes A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}).
#'                        If \code{TRUE}, values \eqn{x_i} satisfying that \eqn{x_i < Q_2 - 10 \cdot \mathrm{IQR}\lbrace x_i \rbrace}
#'                        or \eqn{x_i > Q_3 + 10 \cdot \mathrm{IQR}\lbrace x_i \rbrace} are removed.
#'
#' @return A \code{data.table} object. Contains component summary statistics for each
#'         parameter value combination in \code{parameters} as well as each component.
#'         The following statistics are calculated for each combination:
#'         \itemize{
#'          \item \code{min}: The sample minimum
#'          \item \code{q1}: The empirical 1st percentile.
#'          \item \code{q2.5}: The empirical 2.5th percentile.
#'          \item \code{q5}: The empirical 5th percentile.
#'          \item \code{q10}: The empirical 10th percentile.
#'          \item \code{Q1}: The empirical first quartile of the values (25th percentile)
#'          \item \code{median}: The empirical second quartile of the values (50th percentile)
#'          \item \code{Q3}: The empirical third quartile of the values (75th percentile)
#'          \item \code{q0.90}: The empirical 90th percentile.
#'          \item \code{q0.95}: The empirical 95th percentile.
#'          \item \code{q0.975}: The empirical 97.5th percentile.
#'          \item \code{q0.99}: The empirical 99th percentile.
#'          \item \code{max}: The sample maximum
#'          \item \code{mean}: The sample mean
#'          \item \code{var}: The sample variance
#'          \item \code{sd}: The sample standard deviation
#'          \item \code{mad}: The sample median absolute deviation
#'          \item \code{skewness}: The sample skewness
#'          \item \code{kurtosis}: The sample excess kurtosis
#'         }
#' @export
#'
#' @examples print(1)
simulate_components_statistics <- function(parameters,
                                           N = 100,
                                           shift = FALSE,
                                           attempt_fast = FALSE,
                                           include_pe = FALSE,
                                           parallel = TRUE,
                                           max_cores = 4,
                                           seed = 99,
                                           percent = TRUE,
                                           remove_extremes = FALSE){
  value <- outlier <- NULL
  multiplier <- ifelse(percent, 100, 1)
  simulated_components <- replicate_simulate_components(parameters = parameters,
                                                        N = N,
                                                        shift = shift,
                                                        attempt_fast = attempt_fast,
                                                        include_pe = include_pe,
                                                        parallel = parallel,
                                                        max_cores = max_cores,
                                                        seed = seed)
  if(isTRUE(remove_extremes)){
    simulated_components[, outlier := value < quantile(value, probs = 0.25, na.rm = TRUE) - 10 * IQR(value, na.rm = TRUE) | value > quantile(value, probs = 0.75, na.rm = TRUE) + 10 * IQR(value, na.rm = TRUE), by = c(names(parameters), "component")]
    simulated_components_statistics <- simulated_components[outlier == FALSE, list(min = min(value, na.rm = TRUE),
                                                                                   q1 = quantile(value, probs = 0.01, na.rm = TRUE),
                                                                                   q2.5 = quantile(value, probs = 0.025, na.rm = TRUE),
                                                                                   q5 = quantile(value, probs = 0.05, na.rm = TRUE),
                                                                                   q10 = quantile(value, probs = 0.10, na.rm = TRUE),
                                                                                   Q1 = quantile(value, probs = 0.25, na.rm = TRUE),
                                                                                   median = median(value, na.rm = TRUE),
                                                                                   Q3 = quantile(value, probs = 0.75, na.rm = TRUE),
                                                                                   q90 = quantile(value, probs = 0.90, na.rm = TRUE),
                                                                                   q95 = quantile(value, probs = 0.95, na.rm = TRUE),
                                                                                   q97.5 = quantile(value, probs = 0.975, na.rm = TRUE),
                                                                                   q99 = quantile(value, probs = 0.99, na.rm = TRUE),
                                                                                   max = max(value, na.rm = TRUE),
                                                                                   mean = mean(value, na.rm = TRUE),
                                                                                   var = var(value, na.rm = TRUE),
                                                                                   sd = sd(value, na.rm = TRUE),
                                                                                   mad = mad(value, na.rm = TRUE),
                                                                                   skewness = skewness(value),
                                                                                   kurtosis = kurtosis(value)),
                                                            by = c(names(parameters), "component")]
  }
  else{
    simulated_components_statistics <- simulated_components[, list(min = min(value, na.rm = TRUE),
                                                                   q1 = quantile(value, probs = 0.01, na.rm = TRUE),
                                                                   q2.5 = quantile(value, probs = 0.025, na.rm = TRUE),
                                                                   q5 = quantile(value, probs = 0.05, na.rm = TRUE),
                                                                   q10 = quantile(value, probs = 0.10, na.rm = TRUE),
                                                                   Q1 = quantile(value, probs = 0.25, na.rm = TRUE),
                                                                   median = median(value, na.rm = TRUE),
                                                                   Q3 = quantile(value, probs = 0.75, na.rm = TRUE),
                                                                   q90 = quantile(value, probs = 0.90, na.rm = TRUE),
                                                                   q95 = quantile(value, probs = 0.95, na.rm = TRUE),
                                                                   q97.5 = quantile(value, probs = 0.975, na.rm = TRUE),
                                                                   q99 = quantile(value, probs = 0.99, na.rm = TRUE),
                                                                   max = max(value, na.rm = TRUE),
                                                                   mean = mean(value, na.rm = TRUE),
                                                                   var = var(value, na.rm = TRUE),
                                                                   sd = sd(value, na.rm = TRUE),
                                                                   mad = mad(value, na.rm = TRUE),
                                                                   skewness = skewness(value),
                                                                   kurtosis = kurtosis(value)),
                                                            by = c(names(parameters), "component")]
  }

  if(any("cvx" == names(simulated_components_statistics))){
    simulated_components_statistics$cvx <- simulated_components_statistics$cvx * multiplier
  }
  if(any("cvy" == names(simulated_components_statistics))){
    simulated_components_statistics$cvy <- simulated_components_statistics$cvy * multiplier
  }
  if(any("qran" == names(simulated_components_statistics))){
    simulated_components_statistics$qran <- simulated_components_statistics$qran * multiplier
  }
  return(simulated_components_statistics)

}

