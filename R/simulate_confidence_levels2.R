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
#' @param include_deming A \code{logical} value. If \code{TRUE}, the inside rates based on Deming regression are also included.
#' @param include_ci A \code{logical} value. If \code{TRUE}, approximate confidence intervals for are calculated.
#' @param parallel A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, simulations are done on multiple cores on your local machine for the duration of the function call.
#' @param max_cores An \code{integer} value between \code{2} and the number of cores you have on your computer. If you wish to use fewer cores than is on your computer, you should specify this here.
#' @param seed An \code{integer} value. If you wish to set a seed for reproducibility, set it here. If you do not want reproducibility, set this to \code{NULL}.
#' @param percent A non-missing \code{logical} value (\code{TRUE} / \code{FALSE}). If set to \code{TRUE}, all numbers that are possible to present as percentages are presented as percentages.
#'
#' @description
#' Simulates Monte Carlo approximated empirical confidence level, \eqn{\theta}
#' based on every unique combination of parameters found in \code{parameters},
#' based on \code{N} inside checks.
#'
#' @details
#' The approximate smoothing spline prediction interval for a new observation
#' \eqn{y_{n+1}} given the new predictor value \eqn{x_{n+1}} and a fitted smoothed
#' spline \eqn{\hat{f}} is defined to be:
#'
#' \eqn{PI_T(y_{n+1} | x_{n+1}) \approx \hat{f}(x_{n+1}) \pm t_{n - \mathrm{df}(c), 1 - \alpha/2} \sqrt{\hat{\sigma} + \widehat{\mathrm{Var}}(\hat{f}(x_{n+1}))}}.
#'
#' This function repeats the following steps \code{N} times for each unique parameter
#' combination in \code{parameters}:
#' \enumerate{
#'  \item Simulate \eqn{(x_i, y_i)} for \eqn{i = 1, 2, \ldots, n, n+1}
#'  \item Fit a smoothing spline using \eqn{(x_i, y_i)} for \eqn{i = 1, 2, \ldots, n}
#'  \item Calculate \eqn{z_j = \mathbb{I}(y_{n+1} \in PI_T(y_{n+1} | x_{n+1}))}
#' }
#' The empirical confidence level is calculated using \eqn{\hat{\theta} = \frac{\sum_{j=1}^{N}z_j}{N}}.
#' Assuming that \eqn{z_j \sim \mathrm{Bernoulli}(\theta)}, an asymptotically efficient estimator
#' for the variance of \eqn{\hat{\theta}} is
#'
#' \eqn{\widehat{\mathrm{Var}}[\hat{\theta}] = \frac{\hat{\theta}(1-\hat{\theta})}{N}}
#'
#' The corresponding standard error is denoted as \code{sem} in the output of this function.
#' An asymptotic \eqn{95\%} confidence interval for \eqn{\theta} is
#'
#' \eqn{KI_\theta = \hat{\theta} \pm 1.96 \cdot \mathrm{sem}}
#'
#' This asymptotic confidence interval is only calculated if \code{include_ci} is
#' set to \code{TRUE}.
#'
#' @return A \code{data.table}.
#' @export
#'
#' @examples print(1)
simulate_confidence_levels2 <- function(parameters,
                                        N = 100,
                                        m = 1,
                                        level = 0.95,
                                        extrapolate_ok = FALSE,
                                        shift = FALSE,
                                        attempt_fast = FALSE,
                                        include_deming = FALSE,
                                        include_ci = FALSE,
                                        parallel = TRUE,
                                        max_cores = 4,
                                        seed = 99,
                                        percent = TRUE){
  inside <- NULL
  multiplier <- ifelse(percent, 100, 1)
  simulated_insides <- simulate_insides2(parameters = parameters,
                                         N = N,
                                         m = m,
                                         level = level,
                                         extrapolate_ok = extrapolate_ok,
                                         shift = shift,
                                         attempt_fast = attempt_fast,
                                         include_deming = include_deming,
                                         parallel = parallel,
                                         max_cores = max_cores,
                                         seed = seed)
  simulated_confidence_levels <- simulated_insides[, list(confidence_level = mean(inside, na.rm = TRUE) * multiplier,
                                                          number_of_observations = sum(!is.na(inside)),
                                                          sem = sd(inside, na.rm = TRUE) / sqrt(sum(!is.na(inside))) * multiplier,
                                                          lwr = mean(inside, na.rm = TRUE) * multiplier - qnorm(0.975) * sd(inside, na.rm = TRUE) / sqrt(sum(!is.na(inside))) * multiplier,
                                                          upr = min(mean(inside, na.rm = TRUE) * multiplier + qnorm(0.975) * sd(inside, na.rm = TRUE) / sqrt(sum(!is.na(inside))) * multiplier, 1.0 * multiplier)),
                                                   by = c(names(parameters), if(include_deming){"a"}else{NULL})]

  if(!include_ci){
    simulated_confidence_levels$lwr <- simulated_confidence_levels$upr <- NULL
  }

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
