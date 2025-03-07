#' @title Simulate Empirical Confidence Levels
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
#' @param include_ci A non-missing \code{logical} value. If \code{TRUE},
#'                   approximate confidence intervals for \eqn{\theta}
#'                   (theoretical confidence level) are included in the output.
#' @param parallel A non-missing \code{logical} value. If \code{TRUE},
#'                 simulations are done in parallel. May reduce runtime.
#' @param max_cores An \code{integer}. Must be between \code{2} and the number
#'                  of cores you have on your computer. If fewer cores than
#'                  the total number of cores on your computer is desired,
#'                  specify this here.
#' @param seed An \code{integer}. Ensuring reproducibility even if
#'             \code{parallel = TRUE}. If reproducibility is unimportant, set
#'             this to \code{NULL}.
#' @param percent A non-missing \code{logical} value. If \code{TRUE}, all
#'                results possible to present as percentages are presented as
#'                percentages.
#'
#' @description
#' Approximates the theoretical confidence level, \eqn{\theta},
#' using Monte Carlo simulation for every unique parameter combination in
#' \code{parameters}. These approximations are based on \code{N} inside checks
#' for each parameter combination.
#'
#' @details
#' The approximate smoothing spline prediction interval for a new observation
#' \eqn{y_{n+1}} given the new predictor value \eqn{x_{n+1}} and a fitted smoothed
#' spline \eqn{\hat{f}} is:
#'
#' \eqn{PI_T(y_{n+1} | x_{n+1}) \approx \hat{f}(x_{n+1}) \pm
#' t_{n - \mathrm{df}(c), 1 - \alpha/2}
#' \sqrt{\hat{\sigma} + \widehat{\mathrm{Var}}(\hat{f}(x_{n+1}))}}.
#'
#' This function repeats the following steps \code{N} times for each unique
#' parameter combination:
#' \enumerate{
#'  \item Simulate \eqn{(x_i, y_i)} for \eqn{i = 1, 2, \ldots, n, n+1}
#'  \item Fit a smoothing spline using \eqn{(x_i, y_i)} for
#'        \eqn{i = 1, 2, \ldots, n}
#'  \item Calculate \eqn{z_j = \mathbb{I}(y_{n+1} \in PI_T(y_{n+1} | x_{n+1}))}
#' }
#' The empirical confidence level is calculated using
#' \eqn{\hat{\theta} = \frac{\sum_{j=1}^{N}z_j}{N}}.
#' Assuming that \eqn{z_j \sim \mathrm{Bernoulli}(\theta)}, an asymptotically
#' efficient estimator for the variance of \eqn{\hat{\theta}} is
#'
#' \eqn{\widehat{\mathrm{Var}}[\hat{\theta}] =
#' \frac{\hat{\theta}(1-\hat{\theta})}{N}}
#'
#' The corresponding standard error is denoted as \code{sem} in the output of
#' this function. An asymptotic \eqn{95\%} confidence interval for \eqn{\theta}
#' is
#'
#' \eqn{KI_\theta = \hat{\theta} \pm 1.96 \cdot \mathrm{sem}}
#'
#' This asymptotic confidence interval is only calculated if
#' \code{include_ci = TRUE}.
#'
#' General note: Any function named \code{simulate_***}, is generally not
#' intended to be used by end-users. They are used in research and validity
#' testing.
#'
#' @return
#' A \code{data.table}. The approximated confidence level for each unique
#' combination of parameters in \code{parameters}. The resulting
#' \code{data.table} is expected to have
#' \code{N} \eqn{\times} \code{m} \eqn{\times} \code{length(obs_tau)} rows.
#'
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
#' # Approximate confidence levels for smoothing spline prediction intervals
#' # Takes approximately one minute to run
#' output <- simulate_confidence_levels2(parameters = parameters,
#'                                       N = 1e3,
#'                                       m = 1,
#'                                       level = 0.95,
#'                                       parallel = FALSE,
#'                                       percent = TRUE)
#'
#' # Rounding results for cleaner output
#' output$confidence_level <- round(output$confidence_level, 2L)
#' output$sem <- round(output$sem, 2L)
#'
#' # The output
#' print(output)
#'
#'
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
