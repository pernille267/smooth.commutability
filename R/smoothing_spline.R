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
  interior_knots <- quantile(x = unique_x, probs = seq(from = 0, to = 1, length.out = n_interior_knots + 2))[-c(1, n_interior_knots + 2)]
  names(interior_knots) <- NULL
  return(interior_knots)
}

#' Calculate Smooth Weights
#'
#' This function is used internally for calculating smoother weights.
#' @keywords internal
weight_function <- function(data, weight_data, impr_data, lower = 0.25, df = 5, inv = FALSE, output_type = c("vector", "table"), plot_weight_function = FALSE){
  mean_weights <- impr_data$Var_A + impr_data$Var_B
  inv_weights <- weight_data$MP_A + weight_data$MP_B
  inv_weights <- ifelse(is.na(inv_weights), 0, inv_weights)
  inv_weights <- ifelse(inv_weights <= lower * mean_weights, lower * mean_weights, inv_weights)
  mean_mor <- 0.5 * data$MP_A + 0.5 * data$MP_B
  sorted_mean_mor <- sort(mean_mor)
  weight_order <- sapply(mean_mor, function(x) which(x == sorted_mean_mor))
  ss_inv_weights <- smooth.spline(x = mean_mor, y = inv_weights, w = 1 / (mean_mor)^2, df = df, all.knots = TRUE)
  base_weights <- 1 / (inv_weights)
  pred_weights <- predict(ss_inv_weights)$y
  pred_weights <- pred_weights[weight_order]
  if(!inv){
    pred_weights <- 1 / pred_weights
    if(plot_weight_function){
      plot(ss_inv_weights$x, 1/ss_inv_weights$yin)
      lines(ss_inv_weights$x, 1/ss_inv_weights$y, type = "l")
      abline(h = 1/(mean_weights * lower), col = "gray", lty = 2)
    }

  }
  else{
    if(plot_weight_function){
      plot(ss_inv_weights$x, ss_inv_weights$yin)
      lines(ss_inv_weights$x, ss_inv_weights$y, type = "l")
      abline(h = mean_weights * lower, col = "gray", lty = 2)
    }
  }
  if(output_type[1] == "vector"){
    return(pred_weights)
  }
  else{
    weight_data$MP_A <- 0.5 * pred_weights
    weight_data$MP_B <- 0.5 * pred_weights
    return(weight_data)
  }
}

#' Calculate Smoothing Parameters
#'
#' This function is used internally for smoothing parameters of \code{smoothing_spline()}.
#' @keywords internal
get_smoothing_parameters <- function(y, weights, B, BTB, BTy, Omega, df, lambda, df_max){

  # Calculate ratio of matrices
  r <- sum(diag(BTB)[-c(1:2, (nrow(BTB) - 2):nrow(BTB))]) / sum(diag(Omega)[-c(1:2, (nrow(Omega) - 2):nrow(Omega))])

  # If both df and lambda is null => numeric opt.
  if(is.null(df) && is.null(lambda)){

    nfold_cv_opt <- optimise(f = nfold_cv, lower = 0.38, upper = 1.5, y = y, weights = weights, B = B, BTB = BTB, Omega = Omega, cross_validation = TRUE, r = r)
    lambda <- r * 256**(3 * nfold_cv_opt$minimum - 1)
    df0 <- sum(diag(B %*% solve(BTB + lambda * Omega) %*% t(B) %*% diag(weights)))
    if(is.null(df_max)){
      df_max <- length(y)
    }
    if(df0 > df_max){
      lambda <- to_lambda(df = df_max, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
      if(isTRUE(is.na(lambda))){
        stop("Could not find a suitable value of lambda that corresponds with df = ", df_max)
      }
      cv_crit <- nfold_cv(sp = lambda, y = y, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
    }
    else{
      cv_crit <- nfold_cv_opt$objective
    }
  }

  # If df exists and lambda is NULL => use df
  else if(!is.null(df) && is.null(lambda)){
    lambda <- to_lambda(df = df, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
    if(isTRUE(is.na(lambda))){
      stop("Could not find a suitable value of lambda that corresponds with df = ", df)
    }
    cv_crit <- nfold_cv(sp = lambda, y = y, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
  }

  # If lambda exists and df is NULL => use lambda
  else if(is.null(df) && !is.null(lambda)){
    cv_crit <- nfold_cv(sp = lambda, y = y, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
  }

  # If both lambda and df exist => use lambda
  else if(!is.null(df) && !is.null(lambda)){
    cv_crit <- nfold_cv(sp = lambda, y = y, weights = weights, B = B, BTB = BTB, Omega = Omega, r = r)
  }

  # Gather output in a list object.
  output <- list("lambda" = lambda,
                 "cv_crit" = cv_crit,
                 "r" = r)

  return(output)

}

#' Internal Function for Faster Smoothing Spline Fitting
#'
#' This function is used internally for fitting a smoothing spline faster.
#' @keywords internal
fast_smoothing_spline <- function(orig_x, y, x, weights, all_knots, df = NULL, lambda = NULL, df_max = 7.5){

  # Defaults
  over_ride_df_max <- !is.null(df) | !is.null(lambda)
  control_sp <- list(low = 0.38, high = 1.5)
  smooth.spline_model <- NULL
  W <- diag(weights)

  # If both df and lambda is null => numeric opt.
  if(is.null(df) && is.null(lambda)){
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = orig_x, y = y, w = weights, cv = TRUE, all.knots = all_knots, control.spar = control_sp, keep.stuff = TRUE),
                                    error = function(e) NULL)
    if(is.null(smooth.spline_model)){
      return(NULL)
    }
    df0 <- smooth.spline_model$df
    if(is.null(df_max)){
      df_max <- length(smooth.spline_model$n)
    }
    if(df0 > df_max){
      smooth.spline_model <- tryCatch(expr = smooth.spline(x = orig_x, y = y, w = weights, cv = TRUE, df = df_max, all.knots = all_knots, keep.stuff = TRUE),
                                      error = function(e) NULL)
      if(is.null(smooth.spline_model)){
        return(NULL)
      }
    }
  }

  # If df exists and lambda is NULL => use df
  else if(!is.null(df) && is.null(lambda)){
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = orig_x, y = y, w = weights, cv = TRUE, df = df, all.knots = all_knots, keep.stuff = TRUE),
                                    error = function(e) NULL)
    if(is.null(smooth.spline_model)){
      return(NULL)
    }
  }

  # If lambda exists and df is NULL => use lambda
  else if(is.null(df) && !is.null(lambda)){
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = orig_x, y = y, w = weights, cv = TRUE, lambda = lambda, all.knots = all_knots, keep.stuff = TRUE),
                                    error = function(e) NULL)
  }
  # If both lambda and df exist => use lambda
  else{
    smooth.spline_model <- tryCatch(expr = smooth.spline(x = orig_x, y = y, w = weights, cv = TRUE, lambda = lambda, all.knots = all_knots, keep.stuff = TRUE),
                                    error = function(e) NULL)
  }

  # If smooth.spline() returns an error, return NULL
  if(is.null(smooth.spline_model)){
    return(NULL)
  }

  # Smoothing spline fitting components
  internal_matrices <- get_matrices(smooth.spline_model$auxM)
  B <- spline.des(knots = smooth.spline_model$fit$knot, x = x, derivs = rep(0, smooth.spline_model$n), outer.ok = TRUE)$design
  BTB <- internal_matrices$BTB
  BTy <- internal_matrices$BTy
  BTWWB <- t(B) %*% tcrossprod(W, W) %*% B
  Omega <- internal_matrices$Omega
  Q <- solve(BTB + smooth.spline_model$lambda * Omega)
  beta <- smooth.spline_model$fit$coef
  S <- B %*% Q %*% t(B) %*% W
  df <- smooth.spline_model$df
  fitted <- smooth.spline_model$y
  residuals <- smooth.spline_model$yin - fitted
  var_eps <- smooth.spline_model$pen.crit / (smooth.spline_model$n - df)
  cov_beta <- var_eps * Q %*% BTB %*% Q
  #cov_beta <- var_eps * Q %*% BTWWB %*% Q
  var_fit <- diag(B %*% cov_beta %*% t(B))
  #var_fit <- var_eps * diag(tcrossprod(S, S))
  sp <- smooth.spline_model$spar
  lambda <- smooth.spline_model$lambda

  output <- list(fitted = fitted,
                 residuals = residuals,
                 df = df,
                 df_max = if(over_ride_df_max){"IGNORED"}else{df_max},
                 lambda = lambda,
                 sp = sp,
                 cv_crit = smooth.spline_model$cv.crit,
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
                 weights = weights,
                 x = orig_x,
                 y = y)

  return(output)
}


#' Fit a Smoothing Spline Model with O'Sullivan Penalty
#'
#' This function fits a smoothing spline model to the provided data using an O'Sullivan penalty matrix. It is designed to work with various data types and offers flexibility in model tuning through degrees of freedom (`df`) and penalty parameter (`lambda`).
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

smoothing_spline <- function(data, weights = 1, df = NULL, lambda = NULL, df_max = 7.5, attempt_fast = FALSE, na_rm = TRUE){

  MP_B <- NULL
  cv_crit <- NA
  copy_data <- copy(convert_to_data_table(data))
  over_ride_df_max <- !is.null(df) | !is.null(lambda)
  weights_order <- order(data$MP_B)

  # Check for required variables in data
  if(!any("MP_B" == names(copy_data)) | !any("MP_A" == names(copy_data))){
    stop("'data' is expected to have variables named 'MP_B' and 'MP_A'.")
  }

  if(isTRUE(na_rm)){
    if(any(is.na(copy_data$MP_B)) | any(is.na(copy_data$MP_A))){
      copy_data <- na.omit(copy_data)
    }
  }

  # Sort according to predictor values
  setorder(x = copy_data, MP_B)
  orig_x <- copy_data$MP_B
  y <- copy_data$MP_A

  # Transform predictor values to [0, 1] to avoid overflow
  x <- (orig_x - orig_x[1]) / diff(range(orig_x))

  # Calculation of interior knots
  interior_knots <- calculate_interior_knots(x)

  # Calculation and handling of weights
  if(is.null(weights) | length(weights) == 1){
    weights <- rep(1, length(x))
  }
  else{
    weights <- weights[weights_order]
    weights <- weights / mean(weights)
  }

  # Faster fitting of the smoothing spline, using smooth.spline() functionality
  if(isTRUE(attempt_fast)){
    attempted_fast <- fast_smoothing_spline(orig_x = orig_x, y = y, x = x, weights = weights, all_knots = c(0, interior_knots, 1), df = df, lambda = lambda, df_max = df_max)
    if(is.null(attempted_fast)){
      Recall(data = copy_data, weights = weights, df = df, lambda = lambda, df_max = df_max, attempt_fast = FALSE, na_rm = na_rm)
    }
    output <- attempted_fast
    class(output) <- "smoothing_spline"
    return(invisible(output))
  }

  # Make diagonal matrix from weights
  W <- diag(weights)

  # Calculation of boundary knots and B-spline basis matrix
  x_min_max <- range(x)
  B <- spline.des(knots = c(rep(0, 4), interior_knots, rep(1, 4)), x = x, derivs = rep(0, length(x)), outer.ok = TRUE)$design

  # Penalty matrix and cross-products
  Omega <- penalty_matrix(x_min = x_min_max[1], x_max = x_min_max[2], interior_knots = interior_knots)

  # Relevnat matrix products
  BTB <- t(B) %*% W %*% B
  BTy <- t(B) %*% W %*% y
  BTWWB <- t(B) %*% tcrossprod(W, W) %*% B

  # Degrees of freedom / smoothing parameter calculations
  smoothing_pars <- get_smoothing_parameters(y, weights, B, BTB, BTy, Omega, df, lambda, df_max)

  # Extract relevant smoothing parameters
  lambda <- smoothing_pars$lambda
  cv_crit <- smoothing_pars$cv_crit
  r <- smoothing_pars$r

  # Spline fitting
  Q <- solve(BTB + lambda * Omega) # Q is symmetric!
  beta <- Q %*% BTy
  S <- B %*% Q %*% t(B) %*% W
  df <- sum(diag(S))
  fitted <- S %*% y
  residuals <- y - fitted
  var_eps <- sum(weights * residuals^2) / (length(x) - df)
  cov_beta <- var_eps * Q %*% BTB %*% Q
  #cov_beta <- var_eps * Q %*% BTWWB %*% Q
  var_fit <- diag(B %*% cov_beta %*% t(B))
  #var_fit <- var_eps * diag(tcrossprod(S, S))
  sp <- 1/(3 * log(256)) * log(lambda) - 1/(3 * log(256)) * log(r) + 1/3

  output <- list(fitted = fitted,
                 residuals = residuals,
                 df = df,
                 df_max = if(over_ride_df_max){"IGNORED"}else{df_max},
                 lambda = lambda,
                 sp = sp,
                 cv_crit = cv_crit,
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
                 weights = weights,
                 x = orig_x,
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
  cat("LOOCV:", format(output$cv_crit, digits=7))
  cat("\n------------------------------------------------------")
  invisible(output)
}

#' Default Plot Method for Smoothing Spline Output
#'
#' @param output A \code{smoothing_spline} object, which is the output of \code{smoothing_spline()} method.
#' @param error A \code{character} string value. Can be one of the following:
#'              \itemize{
#'                \item \code{confidence}: Include estimated confidence bands around the smoothing spline.
#'                \item \code{prediction}: Include estimated prediction bands around the smoothing spline.
#'                \item \code{none}: Do not include estimated uncertainty bands around the smoothing spline.
#'              }
#' @param level A \code{double}. The confidence level of the error interval.
#'              Must be between \code{0} and \code{1}. Default value is \code{0.95}.
#'
#' @return Returns a \code{ggplot2} object. This plot contains the fitted \code{smoothing_spline} object,
#'         the data points used in fitting it, and uncertainty bands
#'         around the smoothing spline.
#' @export
#' @method plot smoothing_spline
#'
#' @examples
#' print(1)
plot.smoothing_spline <- function(output, error = c("confidence", "prediction", "none"), level = 0.95){
  x <- y <- lwr <- upr <- fitted <- NULL
  fit_variance <- output$var_fit
  pred_variance <- output$var_fit + output$var_eps
  curve_variance <- 0
  plot_main <- "Smoothing Spline Fit:"
  if(error[1] == "confidence"){
    curve_variance <- fit_variance
    plot_main <- paste0("Smoothing Spline Fit + ", round(level * 100, 2L), "% Confidence Intervals:")
  }
  else if(error[1] == "prediction"){
    curve_variance <- pred_variance
    plot_main <- paste0("Smoothing Spline Fit + ", round(level * 100, 2L), "% Prediction Intervals:")
  }
  else{
    curve_variance <- fit_variance / 100
  }
  t_quant <- qt(p = (1 - level)/2, df = length(output$x) - output$df, lower.tail = FALSE)
  plotting_data <- data.table("x" = output$x,
                              "y" = output$y,
                              "fitted" = output$fitted[,],
                              "lwr" = output$fitted[,] - t_quant * sqrt(curve_variance),
                              "upr" = output$fitted[,] + t_quant * sqrt(curve_variance))

  text_placement_x <- output$x[2]
  text_placement_y <- output$y[length(output$y)]

  ggplot(data = plotting_data) +
    geom_ribbon(mapping = aes(x = x, ymin = lwr, ymax = upr), fill = "green", alpha = 0.3, color = "black", outline.type = "full") +
    geom_line(mapping = aes(x = x, y = fitted)) +
    geom_point(mapping = aes(x = x, y = y), shape = 21, fill = "#55CDEC", color = "black") +
    scale_x_continuous(name = "Predictor values", n.breaks = 8) +
    scale_y_continuous(name = "Fitted values", n.breaks = 8) +
    labs(title = plot_main,
         subtitle = paste("Effective degrees of freedom (df) =", round(output$df, 3L), "and LOOCV =", round(output$cv_crit, 3L))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}



