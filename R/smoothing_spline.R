#' Internal Function for Data Conversion
#'
#' This function is used internally for converting different data types to data.table.
#' @keywords internal
convert_to_data_table <- function(data){
  if(is.list(data)){
    setDT(data)
  }
  else if(is.data.frame(data)){
    as.data.table(data)
  }
  else{
    stop("'data' is expected to be a data.table object but is a ", class(data))
  }
}

#' Calculate interior knots
#'
#' This function is used internally for calculating interior knots.
#' @keywords internal
calculate_interior_knots <- function(x) {
  unique_x <- unique(x)
  n_interior_knots <- length(unique_x)
  if(n_interior_knots >= 50){
    n_interior_knots <- .nknots.smspl(n = n_interior_knots)
  }
  interior_knots <- quantile(x = unique_x,
                             probs = seq(from = 0, to = 1, length.out = n_interior_knots + 2))[-c(1, n_interior_knots + 2)]
  names(interior_knots) <- NULL
  return(interior_knots)
}

#' Calculate Smoothing Parameters
#'
#' This function is used internally for smoothing parameters of \code{smoothing_spline()}.
#' @keywords internal
get_smoothing_parameters <- function(y, weights, B, BTB, BTy, Omega, df, lambda, df_max, method = "loocv", sg = c(0.38, 1.50), smudge = 1.2){

  # Calculate ratio of matrices
  r <- calculate_r(BTB, Omega)

  # If both df and lambda is null => numeric opt.
  if(is.null(df) && is.null(lambda)){
    cv_opt <- NULL
    if(method == "gcv"){
      cv_opt <- optimise(f = gcv2,
                         lower = sg[1],
                         upper = sg[2],
                         y = y,
                         weights = weights,
                         B = B,
                         BTB = BTB,
                         Omega = Omega,
                         cross_validation = TRUE,
                         r = r,
                         smudge = smudge)
    }
    else{
      cv_opt <- optimise(f = nfold_cv2,
                         lower = sg[1],
                         upper = sg[2],
                         y = y,
                         weights = weights,
                         B = B,
                         BTB = BTB,
                         Omega = Omega,
                         cross_validation = TRUE,
                         r = r)
    }

    lambda <- r * 256**(3 * cv_opt$minimum - 1)
    df0 <- calculate_df(lambda,
                        weights,
                        B,
                        BTB,
                        Omega)
    if(is.null(df_max)){
      df_max <- length(y)
    }
    if(df0 > df_max){
      sp <- optimise(f = error_df2,
                     lower = -1.5,
                     upper = 1.5,
                     df = df_max,
                     weights = weights,
                     B = B,
                     BTB = BTB,
                     Omega = Omega,
                     r = r)
      if(sp$objective > 1e-2){
        stop("Could not find a suitable value of lambda that corresponds with df = ", df_max, ".")
      }
      lambda <- r * 256 ** (3 * sp$minimum - 1.0)
      if(method == "gcv"){
        cv_crit <- gcv2(sp = lambda,
                        y = y,
                        weights = weights,
                        B = B,
                        BTB = BTB,
                        Omega = Omega,
                        r = r,
                        smudge = smudge)
      }
      else{
        cv_crit <- nfold_cv2(sp = lambda,
                             y = y,
                             weights = weights,
                             B = B,
                             BTB = BTB,
                             Omega = Omega,
                             r = r)
      }
      method <- paste0("restricted_", method)
    }
    else{
      cv_crit <- cv_opt$objective
    }
  }

  # If df exists and lambda is NULL => use df
  else if(!is.null(df) && is.null(lambda)){
    sp <- optimise(f = error_df2,
                   lower = -1.5,
                   upper = 1.5,
                   df = df,
                   weights = weights,
                   B = B,
                   BTB = BTB,
                   Omega = Omega,
                   r = r)
    if(sp$objective > 1e-2){
      stop("Could not find a suitable value of lambda that corresponds with df = ", df, ".")
    }
    lambda <- r * 256 ** (3 * sp$minimum - 1.0)
    if(method == "gcv"){
      cv_crit <- gcv2(sp = lambda,
                      y = y,
                      weights = weights,
                      B = B,
                      BTB = BTB,
                      Omega = Omega,
                      r = r,
                      smudge = smudge)
    }
    else{
      cv_crit <- nfold_cv2(sp = lambda,
                           y = y,
                           weights = weights,
                           B = B,
                           BTB = BTB,
                           Omega = Omega,
                           r = r)
    }
    method <- "df_solve"
  }

  # If lambda exists and df is NULL => use lambda
  else if(is.null(df) && !is.null(lambda)){
    if(method == "gcv"){
      cv_crit <- gcv2(sp = lambda,
                      y = y,
                      weights = weights,
                      B = B,
                      BTB = BTB,
                      Omega = Omega,
                      r = r,
                      smudge = smudge)
    }
    else{
      cv_crit <- nfold_cv2(sp = lambda,
                           y = y,
                           weights = weights,
                           B = B,
                           BTB = BTB,
                           Omega = Omega,
                           r = r)
    }
    method <- "none"
  }
  # If both lambda and df exist => use lambda
  else if(!is.null(df) && !is.null(lambda)){
    if(method == "gcv"){
      cv_crit <- gcv2(sp = lambda,
                      y = y,
                      weights = weights,
                      B = B,
                      BTB = BTB,
                      Omega = Omega,
                      r = r,
                      smudge = smudge)
    }
    else{
      cv_crit <- nfold_cv2(sp = lambda,
                           y = y,
                           weights = weights,
                           B = B,
                           BTB = BTB,
                           Omega = Omega,
                           r = r)
    }
    method <- "none"
  }

  # Gather output in a list object.
  output <- list("lambda" = lambda,
                 "cv_crit" = cv_crit,
                 "r" = r,
                 "method" = method)

  return(output)

}

#' Internal Function for Faster Smoothing Spline Fitting
#'
#' This function is used internally for fitting a smoothing spline faster.
#' @keywords internal
fast_smoothing_spline <- function(x, y, u, weights, all_knots, df = NULL, lambda = NULL, df_max = 7.5, additional_parameters = list("method" = "gcv", "smudge" = 1.2, "c_star_min" = 0.38, "c_star_max" = 1.50, "tol" = 0.5, "window" = 10, "iter" = 1)){

  # Initialize
  method <- "loocv"
  sp_search_interval <- c(0.38, 1.50)
  cv <- TRUE

  if(any("method" == names(additional_parameters))){
    if(additional_parameters$method == "gcv"){
      cv <- FALSE
      method <- "gcv"
    }
  }

  c_star_min <- 0.38
  c_star_max <- 1.50

  if(any("c_star_min" == names(additional_parameters))){
    c_star_min <- additional_parameters$c_star_min
    sp_search_interval[1] <- c_star_min
  }
  if(any("c_star_max" == names(additional_parameters))){
    c_star_max <- additional_parameters$c_star_max
    sp_search_interval[2] <- c_star_max
  }

  # Defaults
  over_ride_df_max <- !is.null(df) | !is.null(lambda)
  control_sp <- list(low = c_star_min, high = c_star_max)
  smooth.spline_model <- NULL
  W <- diag(weights)

  # If both df and lambda is null => numeric opt.
  if(is.null(df) && is.null(lambda)){
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = x,
                                                         y = y,
                                                         w = weights,
                                                         cv = cv,
                                                         all.knots = all_knots,
                                                         control.spar = control_sp,
                                                         tol = 1e-90,
                                                         keep.stuff = TRUE),
                                    error = function(e) NULL)
    if(is.null(smooth.spline_model)){
      return(NULL)
    }
    df0 <- smooth.spline_model$df
    if(is.null(df_max)){
      df_max <- smooth.spline_model$n
    }
    if(df0 > df_max){
      method <- paste0("restricted_", method)
      smooth.spline_model <- tryCatch(expr = smooth.spline(x = x,
                                                           y = y,
                                                           w = weights,
                                                           cv = cv,
                                                           df = df_max,
                                                           all.knots = all_knots,
                                                           tol = 1e-90,
                                                           keep.stuff = TRUE),
                                      error = function(e) NULL)
      if(is.null(smooth.spline_model)){
        return(NULL)
      }
    }
  }

  # If df exists and lambda is NULL => use df
  else if(!is.null(df) && is.null(lambda)){
    method <- "df_solve"
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = x,
                                                         y = y,
                                                         w = weights,
                                                         cv = cv,
                                                         df = df,
                                                         all.knots = all_knots,
                                                         tol = 1e-90,
                                                         keep.stuff = TRUE),
                                    error = function(e) NULL)
    if(is.null(smooth.spline_model)){
      return(NULL)
    }
  }

  # If lambda exists and df is NULL => use lambda
  else if(is.null(df) && !is.null(lambda)){
    method <- "none"
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = x,
                                                         y = y,
                                                         w = weights,
                                                         cv = cv,
                                                         lambda = lambda,
                                                         all.knots = all_knots,
                                                         tol = 1e-90,
                                                         keep.stuff = TRUE),
                                    error = function(e) NULL)
  }
  # If both lambda and df exist => use lambda
  else{
    method <- "none"
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = x,
                                                         y = y,
                                                         w = weights,
                                                         cv = cv,
                                                         lambda = lambda,
                                                         all.knots = all_knots,
                                                         tol = 1e-90,
                                                         keep.stuff = TRUE),
                                    error = function(e) NULL)
  }

  # If smooth.spline() returns an error, return NULL
  if(is.null(smooth.spline_model)){
    return(NULL)
  }

  # Smoothing spline fitting components
  internal_matrices <- get_matrices(smooth.spline_model$auxM)
  B <- spline.des(knots = smooth.spline_model$fit$knot,
                  x = u,
                  derivs = rep(0, smooth.spline_model$n),
                  outer.ok = TRUE)$design
  BTB <- internal_matrices$BTB
  BTy <- internal_matrices$BTy
  BTWWB <- t(B) %*% tcrossprod(W, W) %*% B
  Omega <- internal_matrices$Omega
  Q <- calculate_Q(smooth.spline_model$lambda, BTB, Omega) # Q is symmetric!
  beta <- smooth.spline_model$fit$coef
  S <- calculate_S2(weights, B, Q)
  df <- smooth.spline_model$df
  fitted <- smooth.spline_model$y
  derivatives <- predict(object = smooth.spline_model,
                         x = smooth.spline_model$x,
                         deriv = 1)$y
  residuals <- smooth.spline_model$yin - fitted
  var_eps <- smooth.spline_model$pen.crit / (smooth.spline_model$n - df)
  cov_beta <- var_eps * calculate_cov_beta(BTB, Q)
  var_fit <- calculate_pred_var(B, cov_beta)[,]
  sp <- smooth.spline_model$spar
  lambda <- smooth.spline_model$lambda

  output <- list(fitted = fitted,
                 residuals = residuals,
                 derivatives = derivatives,
                 df = df,
                 df_max = if(over_ride_df_max){"IGNORED"}else{df_max},
                 lambda = lambda,
                 sp = sp,
                 cv_crit = smooth.spline_model$cv.crit,
                 method = method,
                 sp_search_interval = sp_search_interval,
                 var_eps = var_eps,
                 var_fit = var_fit,
                 coefficients = beta,
                 cov_beta = cov_beta,
                 smoothing_matrix = S,
                 penalty_matrix = Omega,
                 B = B,
                 BTB = BTB,
                 interior_knots = all_knots[all_knots != 0 & all_knots != 1],
                 knot_boundary = c(0, 1),
                 n = smooth.spline_model$n,
                 weights = weights,
                 u = u,
                 x = x,
                 y = y)

  return(output)
}

#' Iterative Estimation of Smoothing Spline Weights
#'
#' This function is used internally for estimating weights.
#' @keywords internal
iterative_weights_estimation <- function(x, y, u, interior_knots, df, lambda, df_max, tol = 0.5, window = 10, iter = 1){

  # Initialize weights
  weights0 <- rep(1, length(u))

  # Add small error to 'x' to avoid smooth.spline warning
  x <- x + rnorm(length(x), mean = 0, sd = mean(x, na.rm = TRUE) * 1e-10)

  # Iteratize Algorithm Starts ...

  # STEP 1
  fit0 <- fast_smoothing_spline(x = x,
                                y = y,
                                u = u,
                                weights = weights0,
                                all_knots = c(0, interior_knots, 1),
                                df = df,
                                lambda = lambda,
                                df_max = df_max)

  if(is.null(fit0)){
    weights <- rep(1, length(u))
    return(weights)
  }

  studentized_residuals0 <- fit0$residuals / sqrt(fit0$var_eps * (1 - 1 / length(fit0$u) * fit0$df))
  squared_studentized_residuals0 <- studentized_residuals0 ** 2
  inv_weights1 <- local_average(x = squared_studentized_residuals0,
                                weights = fit0$weights,
                                window = window)

  # New Weights
  weights1 <- 1 / inv_weights1

  # STEP 2
  fit1 <- fast_smoothing_spline(x = x,
                                y = y,
                                u = u,
                                weights = weights1,
                                all_knots = c(0, interior_knots, 1),
                                df = df,
                                lambda = lambda,
                                df_max = df_max)

  if(is.null(fit1)){
    weights <- rep(1, length(u))
    return(weights)
  }

  if(iter <= 1 & fit1$cv_crit / fit0$cv_crit < tol){
    weights <- weights1
    return(weights)
  }
  else if(fit1$cv_crit / fit0$cv_crit >= tol){
    weights <- rep(1, length(u))
    return(weights)
  }

  studentized_residuals1 <- fit1$residuals / sqrt(fit1$var_eps * (1 - 1 / length(fit1$u) * fit1$df))
  squared_studentized_residuals1 <- studentized_residuals1 ** 2
  inv_weights2 <- local_average(x = squared_studentized_residuals1,
                                weights = fit1$weights,
                                window = window)

  # New Weights
  weights2 <- 1 / inv_weights2

  # STEP 3
  fit2 <- fast_smoothing_spline(x = x,
                                y = y,
                                u = u,
                                weights = weights2,
                                all_knots = c(0, interior_knots, 1),
                                df = df,
                                lambda = lambda,
                                df_max = df_max)

  if(is.null(fit2)){
    weights <- weights1
    return(weights)
  }

  if(iter <= 2 & fit2$cv_crit / fit1$cv_crit < tol){
    weights <- weights2
    return(weights)
  }
  else if(fit2$cv_crit / fit1$cv_crit >= tol){
    weights <- weights1
    return(weights)
  }
  # ... Add more iterations is necessary
  weights <- weights2
  return(weights)

}


#' Fit a Smoothing Spline Model with O'Sullivan Penalty
#'
#' @param data A \code{data.table}, \code{data.frame} or \code{list} containing the measurements.
#'             It should include two specific variables named 'MP_B' and 'MP_A'.
#' @param weights A \code{numeric} vector of length equal to the number of unique elements in \code{SampleID}.
#'                Default is \code{1}, which assigns equal weights to every observation.
#' @param df   Optional; a \code{numeric} value specifying the degrees of freedom for the spline model.
#'             Should be between 2 and the number of clinical samples in \code{data}.
#'             The default value (\code{NULL}) triggers automatic selection using n-fold cross-validation.
#' @param lambda Optional; a non-negative \code{numeric} value specifying the penalty strength.
#'               The default (\code{NULL}) results in the function choosing an optimal value based on
#'               n-fold cross-validation to minimize the mean squared prediction error.
#' @param df_max A \code{numeric} value between 2 and the number of clinical samples in \code{data}.
#'               If n-fold cross-validation is used to obtain the optimal effective number of degrees of freedom, what is its upper limit.
#' @param attempt_fast A \code{logical} value. If set to \code{TRUE}, the algorithm will attempt to utilize the \code{smoothing.spline()} method
#'                     to derive \code{df} and \code{lambda}. It is generally not recommended to set this to \code{TRUE}.
#' @param na_rm A \code{logical} value. If set to \code{TRUE}, NA values are removed before fitting the smoothing spline.
#' @param additional_parameters A \code{list} containing additional alternatives for the smoothing spline fit.
#'
#'
#' @description
#' This function fits a smoothing spline model to the provided data using an O'Sullivan penalty matrix.
#' It is designed to work with various data types and offers flexibility in model tuning through degrees
#' of freedom (`df`) and penalty parameter (`lambda`).
#'
#' @return
#' A \code{list} containing the components of the smoothing spline fit, including:
#' \itemize{
#'    \item \code{fitted}: Fitted values of the spline model.
#'    \item \code{residuals}: Residuals from the fitted model.
#'    \item \code{df}: Effective degrees of freedom used.
#'    \item \code{lambda}: Penalty parameter used.
#'    \item \code{sp}: Scale-invariant smoothing-paramater. Equivalent to \code{spar} in smooth.spline
#'    \item \code{cv_crit}: Criterion value from cross-validation.
#'    \item \code{var_eps}: Estimated variance of the model error terms.
#'    \item \code{var_fit}: Estimated variance of the fitted values.
#'    \item \code{coefficients}: Estimated coefficient vector of length \eqn{n + 4}.
#'    \item \code{cov_beta}: Estimated coefficient covariance matrix of dims \eqn{(n + 4) \times (n+4)}.
#'    \item \code{penalty_matrix}: Estimated penalty matrix of dims \eqn{(n + 4) \times (n+4)}.
#'    \item \code{B}: The B-spline matrix evaluated at each \eqn{n} observation and \eqn{n + 4} knots.
#'    \item \code{BTB}: The inner product of the B-spline matrix \code{B}.
#'    \item \code{interior_knots}: The \eqn{n} inner knots (scaled to unit interval)
#'    \item \code{knot_boundary}: The \eqn{2} boundary knots (scaled to unit interval). These will always be 0 and 1.
#'    \item \code{weights}: The used weights for the smoothing spline fit.
#'    \item \code{x}: The original predictor values sorted in increasing order.
#'    \item \code{y}: The original response values sorted in increasing order.
#' }
#'
#'         These components are often utilized in subsequent analyses and visualizations.
#'
#' @export
#'
#' @examples
#' # Example using a dataset with columns MP_B and MP_A
#' library(data.table)
#' data_example <- data.table(MP_B = 1:10, MP_A = rnorm(10))
#' result <- smoothing_spline(data_example, df = 5)
#' print(result$fitted)

smoothing_spline <- function(data, weights = 1, df = NULL, lambda = NULL, df_max = 7.5, attempt_fast = FALSE, na_rm = TRUE, additional_parameters = list("method" = "gcv", "smudge" = 1.2, "c_star_min" = 0.38, "c_star_max" = 1.50, "tol" = 0.5, "window" = 10, "iter" = 1)){

  # Initialization

  # From additional_parameters
  method <- additional_parameters$method
  smudge <- additional_parameters$smudge
  c_star_min <- additional_parameters$c_star_min
  c_star_max <- additional_parameters$c_star_max
  iter <- additional_parameters$iter
  tol <- additional_parameters$tol
  window <- additional_parameters$window

  # Data structuring
  MP_B <- NULL
  cv_crit <- NA
  copy_data <- copy(convert_to_data_table(data))
  over_ride_df_max <- !is.null(df) | !is.null(lambda)
  weights_order <- order(data$MP_B)
  weights_estimated <- FALSE

  # Check for required variables in data
  if(!any("MP_B" == names(copy_data)) | !any("MP_A" == names(copy_data))){
    stop("'data' is expected to have variables named 'MP_B' and 'MP_A'.")
  }

  # Remove NA values if na_rm = TRUE
  if(isTRUE(na_rm)){
    if(any(is.na(copy_data$MP_B)) | any(is.na(copy_data$MP_A))){
      copy_data <- na.omit(copy_data)
    }
  }

  # Sort according to predictor values
  setorder(x = copy_data, MP_B)
  x <- copy_data$MP_B
  y <- copy_data$MP_A

  # Transform predictor values to [0, 1] to avoid overflow
  u <- (x - x[1]) / diff(range(x))

  # Get number of predictor values
  n <- length(u)

  # Calculation of interior knots
  interior_knots <- calculate_interior_knots(u)

  # Calculation and handling of weights
  if(is.character(weights)){
    weights_estimated <- TRUE
    weights <- iterative_weights_estimation(x = x,
                                            y = y,
                                            u = u,
                                            interior_knots = interior_knots,
                                            df = df,
                                            lambda = lambda,
                                            df_max = df_max,
                                            tol = tol,
                                            window = window)
  }
  if(is.null(weights) | length(weights) <= 3){
    weights <- rep(1, n)
  }
  else{
    if(!weights_estimated){
      weights <- weights[weights_order]
    }
    weights <- weights / mean(weights)
  }

  # Faster fitting of the smoothing spline, using smooth.spline() functionality
  if(isTRUE(attempt_fast)){
    attempted_fast <- fast_smoothing_spline(x = x,
                                            y = y,
                                            u = u,
                                            weights = weights,
                                            all_knots = c(0, interior_knots, 1),
                                            df = df,
                                            lambda = lambda,
                                            df_max = df_max,
                                            additional_parameters = additional_parameters)
    if(is.null(attempted_fast)){
      Recall(data = copy_data,
             weights = weights,
             df = df,
             lambda = lambda,
             df_max = df_max,
             attempt_fast = FALSE,
             na_rm = na_rm,
             additional_parameters = additional_parameters)
    }
    output <- attempted_fast
    class(output) <- "smoothing_spline"
    return(invisible(output))
  }

  # Make diagonal matrix from weights
  W <- diag(weights)

  # Calculation of boundary knots and B-spline basis matrix
  x_min_max <- c(0, 1)
  B <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)),
                  x = u,
                  derivs = rep(0, n),
                  outer.ok = TRUE)$design
  B_deriv <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)),
                        x = u,
                        derivs = rep(1, n),
                        outer.ok = TRUE)$design

  # Penalty matrix and cross-products
  Omega <- penalty_matrix(x_min = 0,
                          x_max = 1,
                          interior_knots = interior_knots)

  # Relevevant matrix products
  BTB <- t(B) %*% W %*% B
  BTy <- t(B) %*% W %*% y
  BTWWB <- t(B) %*% tcrossprod(W, W) %*% B

  # Degrees of freedom / smoothing parameter calculations
  smoothing_pars <- get_smoothing_parameters(y = y,
                                             weights = weights,
                                             B = B,
                                             BTB = BTB,
                                             BTy = BTy,
                                             Omega = Omega,
                                             df = df,
                                             lambda = lambda,
                                             df_max = df_max,
                                             method = method,
                                             sg = c(c_star_min, c_star_max),
                                             smudge = smudge)

  # Extract relevant smoothing parameters
  lambda <- smoothing_pars$lambda
  cv_crit <- smoothing_pars$cv_crit
  r <- smoothing_pars$r
  method <- smoothing_pars$method

  # Spline fitting
  Q <- calculate_Q(lambda, BTB, Omega)
  beta <- Q %*% BTy
  S <- calculate_S2(weights, B, Q)
  df <- sum(diag(S))
  fitted <- S %*% y
  derivatives <- B_deriv %*% beta / diff(range(x))
  residuals <- y - fitted
  var_eps <- sum(weights * residuals^2) / (n - df)
  cov_beta <- var_eps * calculate_cov_beta(BTB, Q)
  var_fit <- calculate_pred_var(B, cov_beta)[,]
  sp <- 1/(3 * log(256)) * log(lambda) - 1/(3 * log(256)) * log(r) + 1/3

  output <- list(fitted = fitted,
                 residuals = residuals,
                 derivatives = derivatives,
                 df = df,
                 df_max = if(over_ride_df_max){"IGNORED"}else{df_max},
                 lambda = lambda,
                 sp = sp,
                 cv_crit = cv_crit,
                 method = method,
                 sp_search_interval = c(c_star_min, c_star_max),
                 var_eps = var_eps,
                 var_fit = var_fit,
                 coefficients = beta,
                 cov_beta = cov_beta,
                 smoothing_matrix = S,
                 penalty_matrix = Omega,
                 B = B,
                 BTB = BTB,
                 interior_knots = interior_knots,
                 knot_boundary = x_min_max,
                 n = n,
                 weights = weights,
                 u = u,
                 x = x,
                 y = y)

  class(output) <- "smoothing_spline"

  # Returning the fitted object as a list
  invisible(output)
}

#' Default Print Method for Smoothing Spline Output
#'
#' @param output A \code{smoothing_spline} object, which is the output of \code{smoothing_spline()} method.
#'
#' @return Returns nothing, but prints compact an relevant information about the output of \code{smoothing_spline()}
#' @export
#' @method print smoothing_spline
#'
#' @examples
#' print(1)
print.smoothing_spline <- function(output){

  cat("\nSmoothing Spline Fit:",
      "\n------------------------------------------------------")
  cat("\nScale-invariant Smoothing Parameter sp =", format(output$sp, digits=7),
      "\nComputational lambda =", format(output$lambda, digits=7))
  cat("\n")
  cat("Maximum Degrees of Fredom (df_max):",
      if(is.null(output$df_max) | is.character(output$df_max)){"IGNORED"}else{format(output$df_max, digits=7)},
      "\nEffective Degrees of Freedom (df):", format(output$df, digits=7))
  cat("\n")
  cat("Estimated Residual Variance:", format(output$var_eps, digits=7))
  cat("\n")
  cat("CV Score:", format(output$cv_crit, digits=7))
  cat("\n------------------------------------------------------")
  invisible(output)
}

#' Plot Method for Smoothing Spline Objects
#'
#' @param output An object containing \code{smoothing_spline} fit results with the following components:
#'               \itemize{
#'                \item \code{x:} Predictor values
#'                \item \code{y:} Response values
#'                \item \code{fitted:} Fitted values
#'                \item \code{weights:} Observation weights
#'                \item \code{residuals:} Model residuals
#'                \item \code{var_eps:} Error term variance
#'                \item \code{var_fit:} Fit variance
#'                \item \code{df:} Effective degrees of freedom
#'                \item \code{cv_crit:} Cross-validation score
#'               }
#' @param type \code{Character} string specifying the type of plot. Options are:
#'              \itemize{
#'                \item \code{confidence}: Smoothing spline fit with confidence intervals (default)
#'                \item \code{prediction}: Smoothing spline fit with prediction intervals
#'                \item \code{residual}: Smoothing spline standardized weighted residuals
#'                \item \code{bias}: Bias analysis (required \code{y_true})
#'                \item \code{none}: Basic smoothing spline fit without intervals
#'              }
#' @param level Numeric value between 0 and 1 specifying the confidence/prediction
#'              level (default: \code{0.95}).
#' @param y_true Optional \code{numeric} vector of true response values for bias analysis.
#'
#' @description
#' Creates diagnostic plots for smoothing spline fits, including confidence intervals,
#' prediction intervals, residual analysis, and bias assessment.
#'
#' @details
#' The function produces different types of diagnostic plots:
#' \itemize{
#'   \item Confidence intervals show uncertainty in the mean response
#'   \item Prediction intervals show uncertainty in individual predictions
#'   \item Residual plots display standardized weighted residuals against fitted values
#'   \item Bias plots show systematic deviations when true values are known
#' }
#' The plots include:
#' \itemize{
#'   \item Green shaded regions for interval estimates
#'   \item Black solid line for the fitted values
#'   \item Blue dashed line for linear fit (where applicable)
#'   \item Red dashed line for reference lines (y=x or y=0)
#'   \item Points showing observed data
#' }
#'
#' @return A ggplot2 object containing the requested diagnostic plot
#'
#' @export
#' @method plot smoothing_spline
#'
#' @examples
#' print(1)
plot.smoothing_spline <- function(output, type = c("confidence", "prediction", "residual", "bias", "none"), level = 0.95, y_true = NULL){

  # Check validity of output
  stopifnot(class(output) == "smoothing_spline")

  # Validate level parameter
  stopifnot(
    is.numeric(level),
    length(level) == 1,
    level > 0,
    level < 1
  )

  # Initialize
  x <- y <- lwr <- upr <- fitted <- x_point <- y_point <- y_line <- NULL
  error_term_variance <- output$var_eps
  fit_variance <- output$var_fit
  pred_variance <- output$var_fit  + error_term_variance / output$weights
  weighted_residuals <- sqrt(output$weights) * output$residuals[,]
  n_fit <- sum(output$weights)
  df_fit <- output$df
  cv_fit <- output$cv_crit
  studentized_residuals <- weighted_residuals / sqrt(error_term_variance * (1 - df_fit / n_fit))
  curve_variance <- 0
  plot_main <- "Smoothing Spline Fit:"
  x_name <- "Predictor Values"
  y_name <- "Response Values"
  updated_plot <- NULL

  if(!is.character(type) | is.null(type)){
    type <- "confidence"
  }
  else if(is.character(type) & length(type) >= 1){
    type <- type[1]
    if(!any(type == c("confidence", "prediction", "residual", "bias", "none"))){
      type <- "confidence"
    }
  }
  else if(is.character(type) & length(type) == 0){
    type <- "confidence"
  }
  else if(is.null(type)){
    type <- "confidence"
  }
  else if(is.na(type)){
    type <- "confidence"
  }

  if(type == "confidence"){
    curve_variance <- fit_variance
    x_point <- output$x
    y_point <- output$y
    y_line <- output$fitted[,]
    plot_main <- paste0("Smoothing Spline Fit + ", round(level * 100, 2L), "% Confidence Intervals:")
  }
  else if(type == "prediction"){
    curve_variance <- pred_variance
    x_point <- output$x
    y_point <- output$y
    y_line <- output$fitted[,]
    plot_main <- paste0("Smoothing Spline Fit + ", round(level * 100, 2L), "% Prediction Intervals:")
  }
  else if(type == "residual"){
    curve_variance <- sqrt(.Machine$double.eps)
    ss_mod <- smooth.spline(x = output$x,
                            y = studentized_residuals, df = 5)
    x_point <- ss_mod$x
    y_point <- ss_mod$yin
    y_line <- ss_mod$y
    plot_main <- "Smoothing Spline Residuals:"
    y_name <- "Standardized Weighted Residuals"
  }
  else if(type == "bias"){
    if(!is.null(y_true) & is.numeric(y_true) & (length(output$y) == length(y_true))){
      curve_variance <- sqrt(.Machine$double.eps)
      point_x <- output$x
      point_y <- (y_true - output$fitted[,])
      line_y <- y_true - output$fitted[,]
      plot_main <- "Smoothing Spline Bias:"
      y_name <- "Bias"
    }
    else{
      curve_variance <- sqrt(.Machine$double.eps)
      ss_mod <- smooth.spline(x = output$x,
                              y = studentized_residuals, df = 5)
      x_point <- ss_mod$x
      y_point <- ss_mod$yin
      y_line <- ss_mod$y
      plot_main <- "Smoothing Spline Residuals:"
      y_name <- "Standardized Weighted Residuals"
    }
  }
  else{
    curve_variance <- sqrt(.Machine$double.eps)
    x_point <- output$x
    y_point <- output$y
    y_line <- output$fitted
  }
  t_quant <- qt(p = (1 - level) / 2.0, df = n_fit - df_fit, lower.tail = FALSE)
  plotting_data <- data.table("x" = x_point,
                              "y" = y_point,
                              "fitted" = y_line,
                              "lwr" = y_line - t_quant * sqrt(curve_variance),
                              "upr" = y_line + t_quant * sqrt(curve_variance))


  base_plot <- ggplot() +
    geom_ribbon(data = plotting_data,
                mapping = aes(x = x, ymin = lwr, ymax = upr),
                fill = "green",
                alpha = 0.3,
                color = "black",
                outline.type = "full") +
    geom_line(data = plotting_data,
              mapping = aes(x = x, y = fitted)) +
    geom_point(data = plotting_data,
               mapping = aes(x = x, y = y),
               shape = 21,
               fill = "#55CDEC",
               color = "black") +
    scale_x_continuous(name = x_name,
                       n.breaks = 8) +
    scale_y_continuous(name = y_name,
                       n.breaks = 8) +
    labs(title = plot_main,
         subtitle = paste("Effective degrees of freedom (df) =",
                          round(df_fit, 3L),
                          "and CV Score =",
                          round(cv_fit, 3L)))

  if(type == "confidence"){
    updated_plot <- base_plot +
      geom_smooth(mapping = aes(x = x, y = y),
                  method = "lm",
                  formula = y ~ x,
                  se = FALSE,
                  color = "blue",
                  linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
  }
  else if(type == "prediction"){
    updated_plot <- base_plot +
      geom_smooth(mapping = aes(x = x, y = y),
                  method = "lm",
                  formula = y ~ x,
                  se = FALSE,
                  color = "blue",
                  linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
  }
  else if(type == "residual"){
    updated_plot <- base_plot +
      geom_abline(slope = 0, intercept = 0, color = "red", linetype = "dashed")
  }
  else if(type == "bias"){
    updated_plot <- base_plot +
      geom_abline(slope = 0, intercept = 0, color = "red", linetype = "dashed")
  }
  else{
    updated_plot <- base_plot +
      geom_smooth(mapping = aes(x = x, y = y),
                  method = "lm",
                  formula = y ~ x,
                  se = FALSE,
                  color = "blue",
                  linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
  }

  if(is.null(updated_plot)){
    stop("Could not render plot...")
  }

  output_plot <- updated_plot +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#000000"),
          plot.subtitle = element_text(hjust = 0.5, color = "#000000"),
          axis.text = element_text(color = "#000000"),
          axis.title = element_text(color = "#000000"))

  output_plot


}



