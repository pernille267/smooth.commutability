#' Simulate empirical confidence levels for a set of parameters
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
#' @param percent A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, all numbers that are possible to present as percentages are presented as percentages.
#'
#' @return A \code{data.table}.
#' @export
#'
#' @examples print(1)
simulate_confidence_levels2 <- function(parameters, N = 100, m = 1, level = 0.99, extrapolate_ok = FALSE, shift = FALSE, attempt_fast = FALSE, parallel = TRUE, max_cores = 4, percent = TRUE){
  multiplier <- ifelse(percent, 100, 1)
  simulated_insides <- simulate_insides2(parameters = parameters,
                                         N = N,
                                         m = m,
                                         level = level,
                                         extrapolate_ok = extrapolate_ok,
                                         shift = shift,
                                         attempt_fast = attempt_fast,
                                         parallel = parallel,
                                         max_cores = max_cores)
  simulated_confidence_levels <- simulated_insides[, list(confidence_level = mean(inside, na.rm = TRUE) * multiplier,
                                                          number_of_observations = sum(!is.na(inside)),
                                                          sem = sd(inside, na.rm = TRUE) / sqrt(sum(!is.na(inside))) * multiplier),
                                                   by = names(parameters)]
  if(any("cvx" == names(simulated_confidence_levels))){
    simulated_confidence_levels$cvx <- simulated_confidence_levels$cvx * multiplier
  }
  if(any("cvy" == names(simulated_confidence_levels))){
    simulated_confidence_levels$cvy <- simulated_confidence_levels$cvy * multiplier
  }
  if(any("qran" == names(simulated_confidence_levels))){
    simulated_confidence_levels$qran <- simulated_confidence_levels$qran * multiplier
  }

  return(simulated_confidence_levels)
}
