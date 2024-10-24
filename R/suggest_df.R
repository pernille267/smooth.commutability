#' Internal Function for Plotting Smoothing Spline Residual Histograms
#'
#' This function is used internally within \code{smoothing_spline_residual_plots()}.
#' @keywords internal
smoothing_spline_residual_histogram <- function(plot_data, studentized = TRUE, df_t = 20){
  stu_res <- raw_res <- NULL
  if(studentized){
    theoretical_density <- data.table("stu_res" = seq(from = min(plot_data$stu_res), to = max(plot_data$stu_res), length.out = 1e2L))
    theoretical_density[, density := dt(x = stu_res, df = df_t)]
    output <- ggplot() +
      geom_histogram(data = plot_data, mapping = aes(x = stu_res, y = after_stat(density)), fill = "blue", color = "black", alpha = 0.5, bins = nclass.Sturges(plot_data$stu_res)) +
      geom_line(data = theoretical_density, mapping = aes(x = stu_res, y = density), color = "red") +
      scale_x_continuous(name = "Studentized residuals", n.breaks = 10) +
      scale_y_continuous(name = "Density", n.breaks = 10) +
      theme_bw()
    return(output)
  }
  output <- ggplot(data = plot_data) +
    geom_histogram(mapping = aes(x = raw_res, y = after_stat(density)), fill = "blue", color = "black", alpha = 0.5, bins = nclass.Sturges(plot_data$raw_res)) +
    scale_x_continuous(name = "Raw residuals", n.breaks = 10) +
    scale_y_continuous(name = "Density", n.breaks = 10) +
    theme_bw()
  return(output)
}

#' Internal Function for Plotting Smoothing Spline Residual Plots
#'
#' This function is used internally within \code{smoothing_spline_residual_plots()}.
#' @keywords internal
smoothing_spline_residual_scatter_plot <- function(plot_data, studentized = TRUE, df_t = 20){
  MP_B <- stu_res <- NULL
  if(studentized){
    one_sd <- pt(q = (-1) * df_t / (df_t - 2), df_t)
    two_sd <- pt(q = (-2) * df_t / (df_t - 2), df_t)
    three_sd <- pt(q = (-3) * df_t / (df_t - 2), df_t)
    output <- ggplot(data = plot_data) +
      geom_ribbon(mapping = aes(x = MP_B, ymin = qt(p = one_sd, df = df_t), ymax = 1 - qt(p = one_sd, df = df_t)), fill = "green", alpha = 0.3, color = "black") +
      geom_ribbon(mapping = aes(x = MP_B, ymin = qt(p = two_sd, df = df_t), ymax = 1 - qt(p = two_sd, df = df_t)), fill = "green", alpha = 0.15, color = "black") +
      geom_ribbon(mapping = aes(x = MP_B, ymin = qt(p = three_sd, df = df_t), ymax = 1 - qt(p = three_sd, df = df_t)), fill = "green", alpha = 0.05, color = "black") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(mapping = aes(x = MP_B, y = stu_res), shape = 21, fill = "blue", alpha = 0.5, color = "black") +
      scale_x_continuous(name = "Predictor values", n.breaks = 10) +
      scale_y_continuous(name = "Studentized residuals", n.breaks = 10) +
      theme_bw()
    return(output)
  }
  output <- ggplot(data = plot_data) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(mapping = aes(x = MP_B, y = stu_res), shape = 21, fill = "blue", alpha = 0.5, color = "black") +
    scale_x_continuous(name = "Predictor values", n.breaks = 10) +
    scale_y_continuous(name = "Raw residuals", n.breaks = 10) +
    theme_bw()
  return(output)
}


#' Plot Smoothing Spline Residual Based plots
#'
#' @param smoothing_spline_object A \code{object} of either class \code{smooth.spline} or \code{smoothing_spline}.
#' @param which_plot A \code{character} string determining which type of residual plot to be plotted. Possible entries include \code{all} (for all four plot variants),
#'                   \code{histogram_studentized}, \code{scatter_studentized}, \code{histogram_raw} and \code{scatter_raw}.
#' @description This function is used to plot residual plots and residual histograms, meant for model diagnosis and validation purposes.
#' @details Note that choosing either \code{scatter_studentized} or \code{histogram_studentized} will use studentized residuals where the
#'          where residuals are scaled by (sig^2 * (1-s_ii))^(-1/2), where s_ii is the leverage for the i-th smallest predictor value.
#'
#' @return A \code{list} of \code{ggplot2} objects if \code{which_plot = "all"}. Otherwise a single \code{ggplot2} object, containing necessary information on the plot.
#' @export
#'
#' @examples
#' # Example using a fictive dataset
#' sim_pars <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10, cve = 0)
#' sim_data <- simulate_eqa_data2(sim_pars, type = 2, AR = TRUE)
#' smoothing_spline_fit <- smooth.spline(x = sim_data$MP_B, y = sim_data$MP_A, df = 5)
#' smoothing_spline_residual_plots(smoothing_spline_fit, "histogram_studentized")
#'
smoothing_spline_residual_plots <- function(smoothing_spline_object, which_plot = c("all", "histogram_studentized", "scatter_studentized", "histogram_raw", "scatter_raw")){
  is_ss_obj_stats <- class(smoothing_spline_object) == "smooth.spline"
  is_ss_obj <- class(smoothing_spline_object) == "smoothing_spline"
  plot_data <- NULL
  df_t <- NULL
  if(!is_ss_obj_stats & !is_ss_obj){
    stop("smoothing_spline_object must either be a smooth.spline or smoothing_spline object. It is neither right now.")
  }
  else if(is_ss_obj_stats){
    df_t <- smoothing_spline_object$n - smoothing_spline_object$df
    var_eps <- smoothing_spline_object$pen.crit / df_t
    raw_res <- smoothing_spline_object$yin - smoothing_spline_object$y
    stu_res <- raw_res / sqrt(var_eps * (1 - smoothing_spline_object$lev))
    plot_data <- data.table("MP_B" = smoothing_spline_object$x,
                            "MP_A" = smoothing_spline_object$yin,
                            "raw_res" = raw_res,
                            "stu_res" = stu_res)
  }
  else if(is_ss_obj){
    df_t <- length(smoothing_spline_object$x) - smoothing_spline_object$df
    var_eps <- smoothing_spline_object$var_eps
    raw_res <- smoothing_spline_object$residuals
    stu_res <- raw_res / sqrt(var_eps * (1 - diag(smoothing_spline_object$S)))
    plot_data <- data.table("MP_B" = smoothing_spline_object$x,
                            "MP_A" = smoothing_spline_object$y,
                            "raw_res" = raw_res,
                            "stu_res" = stu_res)
  }
  histogram_studentized <- NULL
  histogram_raw <- NULL
  scatter_studentized <- NULL
  scatter_raw <- NULL

  if(length(which_plot) > 1){
    which_plot <- which_plot[1]
  }

  if(which_plot == "all"){
    histogram_studentized <- smoothing_spline_residual_histogram(plot_data = plot_data, studentized = TRUE, df_t = df_t)
    histogram_raw <- smoothing_spline_residual_histogram(plot_data = plot_data, studentized = FALSE, df_t = df_t)
    scatter_studentized <- smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = TRUE, df_t = df_t)
    scatter_raw <- smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = FALSE, df_t = df_t)
    return(list("histogram_studentized" = histogram_studentized,
                "histogram_raw" = histogram_raw,
                "scatter_studentized" = scatter_studentized,
                "scatter_raw" = scatter_raw))
  }
  else if(which_plot == "histogram_studentized"){
    return(smoothing_spline_residual_histogram(plot_data = plot_data, studentized = TRUE, df_t = df_t))
  }
  else if(which_plot == "histogram_raw"){
    return(smoothing_spline_residual_histogram(plot_data = plot_data, studentized = FALSE, df_t = df_t))
  }
  else if(which_plot == "scatter_studentized"){
    return(smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = TRUE, df_t = df_t))
  }
  else if(which_plot == "scatter_raw"){
    return(smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = FALSE, df_t = df_t))
  }
  else{
    histogram_studentized <- smoothing_spline_residual_histogram(plot_data = plot_data, studentized = TRUE, df_t = df_t)
    histogram_raw <- smoothing_spline_residual_histogram(plot_data = plot_data, studentized = FALSE, df_t = df_t)
    scatter_studentized <- smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = TRUE, df_t = df_t)
    scatter_raw <- smoothing_spline_residual_scatter_plot(plot_data = plot_data, studentized = FALSE, df_t = df_t)
    return(list("histogram_studentized" = histogram_studentized,
                "histogram_raw" = histogram_raw,
                "scatter_studentized" = scatter_studentized,
                "scatter_raw" = scatter_raw))
  }
}

#' Internal Function for Calculting Zeta Loss With Respect To Effective DOF
#'
#' This function is used internally within multiple functions.
#' @keywords internal
zeta_loss <- function(df, x, y, weights, var_v, var_h, all_knots, loss = "abs", R = 3, mor = TRUE, target = 1){
  if(mor){
    var_v <- var_v / R
    var_h <- var_h / R
  }

  if(length(df) > 1){
    zeta_diffs <- sapply(X = df, FUN = function(df_candidate){
      ss_candidate <- smooth.spline(x = x, y = y, w = weights, df = df_candidate, all.knots = all_knots, cv = FALSE)
      ss_candidate_mean_squared_slope <- mean(predict(ss_candidate, x = x, deriv = 1)$y^2, na.rm = TRUE)
      var_eps <- ss_candidate$pen.crit / (ss_candidate$n - ss_candidate$df)
      if(loss == "sqe"){
        return((var_eps / (var_v + var_h * ss_candidate_mean_squared_slope) - target)^2)
      }
      else{
        return(abs(var_eps / (var_v + var_h * ss_candidate_mean_squared_slope) - target))
      }

    })
    return(zeta_diffs)
  }
  ss_candidate <- smooth.spline(x = x, y = y, w = weights, df = df, all.knots = all_knots, cv = FALSE)
  ss_candidate_mean_squared_slope <- mean(predict(ss_candidate, x = x, deriv = 1)$y^2, na.rm = TRUE)
  var_eps <- ss_candidate$pen.crit / (ss_candidate$n - ss_candidate$df)
  if(loss == "sqe"){
    zeta_diffs <- (var_eps / (var_v + var_h * ss_candidate_mean_squared_slope) - target)^2
  }
  else{
    zeta_diffs <- abs(var_eps / (var_v + var_h * ss_candidate_mean_squared_slope) - target)
  }
  return(zeta_diffs)
}

#' Internal Function for Calculting Zeta Derivative Loss With Respect To Effective DOF
#'
#' This function is used internally within multiple functions.
#' @keywords internal
zeta_loss_deriv <- function(df, x, y, weights, var_v, var_h, all_knots, loss = "abs", R = 3, mor = TRUE, target = 0){
  h <- .Machine$double.eps^(1/(1+2))
  if(mor){
    var_v <- var_v / R
    var_h <- var_h / R
  }

  if(length(df) > 1){
    zeta_derivs <- sapply(df, FUN = function(df_candidate){
      ss_candidate <- smooth.spline(x = x, y = y, w = weights, df = df_candidate, all.knots = all_knots, cv = FALSE)
      ss_candidate_minus <- smooth.spline(x = x, y = y, w = weights, df = df_candidate - h, all.knots = all_knots, cv = FALSE)
      ss_candidate_plus <- smooth.spline(x = x, y = y, w = weights, df = df_candidate + h, all.knots = all_knots, cv = FALSE)

      ss_candidate_mean_squared_slope <- mean(predict(ss_candidate, x = x, deriv = 1)$y^2, na.rm = TRUE)
      ss_candidate_mean_squared_slope_minus <- mean(predict(ss_candidate_minus, x = x, deriv = 1)$y^2, na.rm = TRUE)
      ss_candidate_mean_squared_slope_plus <- mean(predict(ss_candidate_plus, x = x, deriv = 1)$y^2, na.rm = TRUE)

      var_eps <- ss_candidate$pen.crit / (ss_candidate$n - ss_candidate$df)
      var_eps_minus <- ss_candidate_minus$pen.crit / (ss_candidate_minus$n - ss_candidate_minus$df)
      var_eps_plus <- ss_candidate_plus$pen.crit / (ss_candidate_plus$n - ss_candidate_plus$df)

      zeta <- var_eps / (var_v + var_h * ss_candidate_mean_squared_slope)
      zeta_minus <- var_eps_minus / (var_v + var_h * ss_candidate_mean_squared_slope_minus)
      zeta_plus <- var_eps_plus / (var_v + var_h * ss_candidate_mean_squared_slope_plus)

      first_deriv <- (zeta_plus - zeta_minus) / (2 * h)
      second_deriv <- (zeta_plus - 2 * zeta + zeta_minus) / (h ^ 2)

      if(loss == "sqe"){
        return((first_deriv - target)^2)
      }
      else if(loss == "diff"){
        return(first_deriv - target)
      }
      else if(loss == "second_deriv"){
        return(second_deriv)
      }
      else{
        return(abs(first_deriv - target))
      }
    })
    return(zeta_derivs)
  }

  ss_candidate <- smooth.spline(x = x, y = y, w = weights, df = df, all.knots = all_knots, cv = FALSE)
  ss_candidate_minus <- smooth.spline(x = x, y = y, w = weights,  df = df - h, all.knots = all_knots, cv = FALSE)
  ss_candidate_plus <- smooth.spline(x = x, y = y, w = weights, df = df + h, all.knots = all_knots, cv = FALSE)

  ss_candidate_mean_squared_slope <- mean(predict(ss_candidate, x = x, deriv = 1)$y^2, na.rm = TRUE)
  ss_candidate_mean_squared_slope_minus <- mean(predict(ss_candidate_minus, x = x, deriv = 1)$y^2, na.rm = TRUE)
  ss_candidate_mean_squared_slope_plus <- mean(predict(ss_candidate_plus, x = x, deriv = 1)$y^2, na.rm = TRUE)

  var_eps <- ss_candidate$pen.crit / (ss_candidate$n - ss_candidate$df)
  var_eps_minus <- ss_candidate_minus$pen.crit / (ss_candidate_minus$n - ss_candidate_minus$df)
  var_eps_plus <- ss_candidate_plus$pen.crit / (ss_candidate_plus$n - ss_candidate_plus$df)

  zeta <- var_eps / (var_v + var_h * ss_candidate_mean_squared_slope)
  zeta_minus <- var_eps_minus / (var_v + var_h * ss_candidate_mean_squared_slope_minus)
  zeta_plus <- var_eps_plus / (var_v + var_h * ss_candidate_mean_squared_slope_plus)

  first_deriv <- (zeta_plus - zeta_minus) / (2 * h)
  second_deriv <- (zeta_plus - 2 * zeta + zeta_minus) / (h ^ 2)

  if(loss == "sqe"){
    return((first_deriv - target)^2)
  }
  else if(loss == "diff"){
    return(first_deriv - target)
  }
  else if(loss == "second_deriv"){
    return(second_deriv)
  }
  else{
    return(abs(first_deriv - target))
  }
}

#' Internal Function for Plotting Zeta Loss With Respect To Effective DOF
#'
#' This function is used internally within \code{zeta_df_plots()} method.
#' @keywords internal
zeta_df_plot <- function(plot_data, which_plot = c("zeta", "zeta_loss"), df_unit = FALSE, area = NULL){
  plot_data_copy <- copy(plot_data)
  if(length(which_plot) > 1){
    which_plot <- which_plot[1]
  }
  if(df_unit){
    plot_data_copy$df <- (plot_data_copy$df - min(plot_data_copy$df)) / diff(range(plot_data_copy$df))
  }
  if(which_plot == "zeta"){
    if(!is.null(area) & !df_unit){
      plot_data_copy_2_3 <- copy(plot_data_copy)
      plot_data_copy_2_3 <- plot_data_copy_2_3[df >= 2 & df <= 3, ]
      output <- ggplot() +
        geom_hline(yintercept = 1) +
        geom_ribbon(data = plot_data_copy_2_3, mapping = aes(x = df, ymin = 0, ymax = zeta), fill = "green", alpha = 0.5, color = "black", outline.type = "full") +
        geom_line(data = plot_data_copy, mapping = aes(x = df, y = zeta), col = "black", size = 1) +
        geom_label(mapping = aes(x = 2.5, y = mean(plot_data_copy_2_3$zeta) / 2, label = paste(round(area, 3L))), size = 1) +
        scale_x_continuous(n.breaks = 10) +
        scale_y_continuous(n.breaks = 10) +
        theme_bw()
      return(output)
    }
    else{
      output <- ggplot(data = plot_data_copy) +
        geom_hline(yintercept = 1) +
        geom_line(mapping = aes(x = df, y = zeta), col = "black", size = 1) +
        scale_x_continuous(n.breaks = 10) +
        scale_y_continuous(n.breaks = 10) +
        theme_bw()
      return(output)
    }

  }
  else if(which_plot == "zeta_loss"){
    output <- ggplot(data = plot_data_copy) +
      geom_hline(yintercept = 0) +
      geom_line(mapping = aes(x = df, y = zeta_loss), col = "black", size = 1) +
      scale_x_continuous(n.breaks = 10) +
      scale_y_continuous(n.breaks = 10) +
      theme_bw()
    return(output)
  }
}

#' Internal Function for Plotting Zeta Derivative Loss With Respect To Effective DOF
#'
#' This function is used internally within \code{zeta_df_plots()} method.
#' @keywords internal
zeta_deriv_df_plot <- function(plot_data, which_plot = c("angle", "zeta_deriv_loss"), df_unit = FALSE){
  plot_data_copy <- copy(plot_data)
  if(length(which_plot) > 1){
    which_plot <- which_plot[1]
  }
  if(df_unit){
    plot_data_copy$df <- (plot_data_copy$df - min(plot_data_copy$df)) / diff(range(plot_data_copy$df))
  }
  if(which_plot == "angle"){
    output <- ggplot(data = plot_data_copy) +
      geom_hline(yintercept = 0) +
      geom_line(mapping = aes(x = df, y = zeta_raw_angle), col = "black", size = 1) +
      scale_x_continuous(n.breaks = 10) +
      scale_y_continuous(name = "Derivative zeta angle", n.breaks = 10) +
      theme_bw()
    return(output)
  }
  else if(which_plot == "zeta_deriv_loss"){
    output <- ggplot(data = plot_data_copy) +
      geom_hline(yintercept = 0) +
      geom_line(mapping = aes(x = df, y = zeta_loss_deriv), col = "black", size = 1) +
      scale_x_continuous(n.breaks = 10) +
      scale_y_continuous(name = "Loss derivative zeta", n.breaks = 10) +
      theme_bw()
    return(output)
  }
}

#' Internal Function for Plotting Zeta Derivative Loss With Respect To Effective DOF
#'
#' This function is used internally within \code{suggest_df()} method.
#' @keywords internal
zeta_df_plots <- function(df = NULL, x, y, weights, var_v, var_h, all_knots, loss = "abs", R = 3, mor = TRUE, targets = c(1, 0), which_plot = c("all", "zeta", "zeta_loss", "angle", "zeta_deriv_loss"), df_unit = FALSE){
  if(is.null(df)){
    if(mor){
      df_min <- 2
      df_max <- length(x) - length(x) * 0.4
    }
    else{
      df_min <- 2
      df_max <- length(x) - length(x) * (R-1) / R
    }

  }
  else{
    df_min <- df[1]
    df_max <- df[2]
  }

  area_2_3 <- NULL

  if(abs(df_min - 2) < 0.05 & df_max >= 3){
    area_2_3 <- integrate(zeta_loss, lower = 2, upper = 3, x = x, y = y, weights = weights, var_v = var_v, var_h = var_h, all_knots = all_knots, loss = "abs", R = R, mor = mor, target = 0)$value
  }
  df_grid <- seq(from = df_min, to = df_max, length.out = 2e2L)
  zeta_loss_grid <- zeta_loss(df = df_grid, x = x, y = y, weights = weights, var_v = var_v, var_h = var_h, all_knots = all_knots, loss = "abs", R = R, mor = mor, target = targets[1])
  zeta_grid <- zeta_loss(df = df_grid, x = x, y = y, weights = weights, var_v = var_v, var_h = var_h, all_knots = all_knots, loss = "abs", R = R, mor = mor, target = 0)
  df_grid[1] <- df_grid[1] + 0.05
  zeta_loss_deriv_grid <- zeta_loss_deriv(df = df_grid, x = x, y = y, weights = weights, var_v = var_v, var_h = var_h, all_knots = all_knots, loss = "abs", R = R, mor = mor, target = targets[2])
  zeta_loss_angle_grid <- atan(zeta_loss_deriv(df = df_grid, x = x, y = y, weights = weights, var_v = var_v, var_h = var_h, all_knots = all_knots, loss = "diff", R = R, mor = mor, target = targets[2])) * 180 / pi
  df_grid[1] <- df_grid[1] - 0.05
  plot_data <- data.table("df" = df_grid,
                          "zeta" = zeta_grid,
                          "zeta_loss" = zeta_loss_grid,
                          "zeta_raw_angle" = zeta_loss_angle_grid,
                          "zeta_loss_deriv" = zeta_loss_deriv_grid)
  if(length(which_plot) > 1){
    which_plot <- which_plot[1]
  }

  df_versus_zeta <- NULL
  df_versus_zeta_loss <- NULL
  df_versus_zeta_angle <- NULL
  df_versus_zeta_deriv_loss <- NULL

  if(which_plot == "all"){
    df_versus_zeta <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta", df_unit = df_unit, area = area_2_3)
    df_versus_zeta_loss <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta_loss", df_unit = df_unit, area = area_2_3)
    df_versus_zeta_angle <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "angle", df_unit = df_unit)
    df_versus_zeta_deriv_loss <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "zeta_deriv_loss", df_unit = df_unit)
    output <- list("df_versus_zeta" = df_versus_zeta,
                   "df_versus_zeta_loss" = df_versus_zeta_loss,
                   "df_versus_zeta_angle" = df_versus_zeta_angle,
                   "df_versus_zeta_deriv_loss" = df_versus_zeta_deriv_loss)
  }
  else if(which_plot == "zeta"){
    df_versus_zeta <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta", df_unit = df_unit, area = area_2_3)
    output <- df_versus_zeta
  }
  else if(which_plot == "zeta_loss"){
    df_versus_zeta_loss <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta_loss", df_unit = df_unit, area = area_2_3)
    output <- df_versus_zeta_loss
  }
  else if(which_plot == "angle"){
    df_versus_zeta_angle <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "angle", df_unit = df_unit)
    output <- df_versus_zeta_angle
  }
  else if(which_plot == "zeta_deriv_loss"){
    df_versus_zeta_deriv_loss <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "zeta_deriv_loss", df_unit = df_unit)
    output <- df_versus_zeta_deriv_loss
  }
  else{
    df_versus_zeta <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta", df_unit = df_unit, area = area_2_3)
    df_versus_zeta_loss <- zeta_df_plot(plot_data = plot_data, which_plot = "zeta_loss", df_unit = df_unit, area = area_2_3)
    df_versus_zeta_angle <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "angle", df_unit = df_unit)
    df_versus_zeta_deriv_loss <- zeta_deriv_df_plot(plot_data = plot_data, which_plot = "zeta_deriv_loss", df_unit = df_unit)
    output <- list("df_versus_zeta" = df_versus_zeta,
                   "df_versus_zeta_loss" = df_versus_zeta_loss,
                   "df_versus_zeta_angle" = df_versus_zeta_angle,
                   "df_versus_zeta_deriv_loss" = df_versus_zeta_deriv_loss)
  }
  return(output)

}


#' Suggest Optimal Starting Values for Effective Degrees of Freedom
#'
#' @param data A \code{list} or \code{data.table} containing columns \code{SampleID},
#'             \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#' @param ns_score A \code{logical}. If \code{TRUE}, calculates non-linearity scores for
#'                 each method used to derive optimal effective DOF values. See details.
#' @param plots A \code{logical}. If \code{TRUE}, generates model assessment and \eqn{\zeta} dynamics
#'              plots.
#' @param models A \code{logical}. If \code{TRUE}, returns the six optimal \code{smooth.spline}
#'               objects in the output.
#' @param use_weights A code{logical}. If \code{TRUE}, uses weights in the algorithms.
#'                    Weights are automatically determined based on \code{data}.
#' @param na_rm A \code{logical}. If \code{TRUE}, removes NA values before modeling. Note
#'              that this affects algorithm results but not measurement variance or replicate
#'              count calculations.
#'
#' @description This function addresses the challenge of determining suitable values for the
#'              degrees of freedom \code{df} parameter in smoothing splines, particularly for
#'              small sample sizes common in commutability evaluations. It provides
#'              comprehensive data analysis and suggests multiple optimal \code{df} values using
#'              six different algorithms.
#'
#' @details
#' Algorithms
#'
#' The function employs six algorithms to determine optimal \code{df} values:
#' \itemize{
#'    \item Four variants of \eqn{\zeta} minimization
#'    \item Leave-One-Out Cross-Validation (LOOCV)
#'    \item Generalized Cross-Validation (GCV)
#' }
#'
#' While none of these methods are entirely robust for typical sample sizes,
#' they provide valuable insights when used in conjunction with the output plots.
#'
#' Non-linearity Scores
#'
#' When \code{ns_score = TRUE}, the function calculates non-linearity scores based on
#' \itemize{
#'    \item For minimizing \eqn{\zeta(\mathrm{df})}: Score = -1 + \eqn{\int_{2}^{3}\zeta(\mathrm{t}) dt}
#'    \item For minimizing \eqn{\zeta'(\mathrm{df})}: Score = \eqn{-0.001 \cdot \sum_{i=1}^{1000} \zeta'(U_i)}
#'    where \eqn{U_i \sim \mathrm{unif}(2, 3)}
#' }
#'
#' Negative or near-zero scores suggest insufficient non-linearity to justify
#' \code{df} values other than 2. Scores are not calculated for LOOCV and GCV algorithms.
#'
#' Output Plots
#'
#' When \code{plots = TRUE}, four lists of plots are generated:
#' \itemize{
#'    \item \code{residual_scatter_plots}: Studentized residual plots for each optimal \code{df}
#'    \item \code{residual_histogram_plots}: Studentized residual histograms for each optimal \code{df}
#'    \item \code{zeta_vs_df_plots_mor}: \eqn{\zeta} vs. \code{df} relationships using mean of replicates (mor) data
#'    \item \code{zeta_vs_df_plots_ar}: \eqn{\zeta} vs. \code{df} relationships using raw \code{data}
#' }
#'
#' These plots complement the tabular output and can guide decision-making when
#' optimal \code{df} values may not be ideal.
#'
#' @return
#' The output structure varies based on \code{ns_score}, \code{plots} and \code{models} parameters:
#' \itemize{
#'    \item If \code{plots = FALSE} and \code{models = FALSE}: A \code{data.table} (2x6 or 3x6 if \code{ns_score = TRUE})
#'    \item Otherwise: A \code{list} containing elements \code{table}, \code{plots}, and/or \code{models}
#' }
#'
#' @export
#'
#' @examples
#' ## Load packages
#' library(smooth.commutability)
#' library(data.table)
#' library(fasteqa)
#'
#' ## Linear data
#' pars <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, b0 = 0, b1 = 1.05)
#' sim_dat <- simulate_eqa_data(pars) |> setDT()
#' suggest_dfs <- suggest_df(data = sim_dat, ns_score = TRUE, plots = FALSE)
#' print(suggest_dfs)
#'
#' ## Non-linear data
#' pars <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0)
#' sim_dat <- simulate_eqa_data2(pars, type = 2, AR = TRUE) |> setDT()
#' suggest_dfs <- suggest_df(data = sim_dat, ns_score = TRUE, plots = FALSE)
#' print(suggest_dfs)

suggest_df <- function(data, ns_score = FALSE, plots = FALSE, models = FALSE, use_weights = FALSE, na_rm = TRUE){

  # Check for required variables in data
  if(!any("MP_B" == names(data)) | !any("MP_A" == names(data)) | !any("SampleID" == names(data)) | !any("ReplicateID" == names(data))){
    stop("'data' is expected to have variables named 'MP_B' and 'MP_A'.")
  }

  # Calculate imprecision estimates
  impr <- global_precision_estimates(data)

  # Remove NA-values if na_rm = TRUE
  if(isTRUE(na_rm)){
    if(any(is.na(data$MP_B)) | any(is.na(data$MP_A))){
      data <- na.omit(data)
    }
  }

  # Calculate average number of replicates used on each sample
  R_i <- count_samplewise_replicates(data, summary = "none")$R_i
  R <- mean(R_i)

  # Calculate mean of replicates (MOR) data
  mor_data <- fun_of_replicates2(data, "mean") |> setDT()

  # Calculate weight data
  if(use_weights){
    vor_data <- fun_of_replicates2(data, "var") |> setDT()
    weight_data_mor <- weight_function(mor_data, vor_data, impr, output_type = "vector")
    weight_data_ar <- sapply(1:length(mor_data$SampleID), function(x) rep(weight_data_mor[x], R_i[x]), simplify = FALSE) |> unlist()
  }
  else{
    vor_data <- fun_of_replicates(data, "var") |> setDT()
    weight_data_mor <- rep(1, length(mor_data$MP_A))
    weight_data_ar <- rep(1, length(data$MP_A))
  }


  axis_var <- if (impr$lambda < 0.5) {
    c("MP_A", "MP_B", "Var_A", "Var_B")
  } else {
    c("MP_B", "MP_A", "Var_B", "Var_A")
  }

  # Extract data
  setorderv(mor_data, axis_var[1], order = 1)
  x_mor <- mor_data[[axis_var[1]]]
  y_mor <- mor_data[[axis_var[2]]]
  x <- data[[axis_var[1]]]
  y <- data[[axis_var[2]]]
  var_h <- impr[[axis_var[3]]]
  var_v <- impr[[axis_var[4]]]

  # Calculate knots to be used (MOR)
  x_unit_mor <- x_mor[order(x_mor)]
  x_unit_mor <- (x_unit_mor - x_unit_mor[1]) / diff(range(x_unit_mor))
  all_knots_mor <- c(0, calculate_interior_knots(x_unit_mor), 1)

  # Calculate knots to be used (AR)
  x_unit_ar <- x[order(x)]
  x_unit_ar <- (x_unit_ar - x_unit_ar[1]) / diff(range(x_unit_ar))
  all_knots_ar <- c(0, calculate_interior_knots(x_unit_ar), 1)

  # Fit two smooth.spline models using LOOCV and GCV
  ss_loo <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, all.knots = all_knots_mor, cv = TRUE, control.spar = list(low = 0.38, high = 2))
  ss_gcv <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, all.knots = all_knots_mor, cv = FALSE, control.spar = list(low = 0.38, high = 2))

  # Extract optimal degrees of freedom based on LOOCV and GCV
  opt_df_loo <- ss_loo$df
  opt_df_gcv <- ss_gcv$df

  # Approximate second derivatives
  candidate_df_mor <- seq(from = 2, to = length(x_mor) * 0.6, length.out = 200)
  candidate_df_ar <- seq(from = 2, to = length(x) / (R - 0.1 * R), length.out = 200)
  zeta_df_mor <- zeta_loss(df = candidate_df_mor, x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 0)
  zeta_df_ar <- zeta_loss(df = candidate_df_ar, x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 0) #WARNS
  second_deriv_object <- data.table(df_mor = candidate_df_mor,
                                    df_ar = candidate_df_ar,
                                    second_deriv_zeta_mor = predict(smooth.spline(x = candidate_df_mor, y = zeta_df_mor, df = 25), deriv = 2)$y,
                                    second_deriv_zeta_ar = predict(smooth.spline(x = candidate_df_ar, y = zeta_df_ar, df = 25), deriv = 2)$y)

  # Get upper search regions for df
  df_max_mor <- obtain_df_max(df = candidate_df_mor, second_deriv = second_deriv_object$second_deriv_zeta_mor)
  df_max_ar <- obtain_df_max(df = candidate_df_ar, second_deriv = second_deriv_object$second_deriv_zeta_ar)
  df_max_mor_global <- min(candidate_df_mor[candidate_df_mor >= candidate_df_mor[which.max(second_deriv_object$second_deriv_zeta_mor)] & abs(second_deriv_object$second_deriv_zeta_mor) == min(abs(second_deriv_object$second_deriv_zeta_mor))], 2 * df_max_mor)
  df_max_ar_global <- min(candidate_df_ar[candidate_df_ar >= candidate_df_ar[which.max(second_deriv_object$second_deriv_zeta_ar)] & abs(second_deriv_object$second_deriv_zeta_ar) == min(abs(second_deriv_object$second_deriv_zeta_ar))], 2 * df_max_ar)

  # Extract optimal degrees of freedom based on constrained and unconstrained zeta loss and zeta deriv loss using MOR data
  opt_df_zm_mor <- optimise(f = zeta_loss, interval = c(2, df_max_mor), x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 1)
  opt_df_zd_mor <- optimise(f = zeta_loss_deriv, interval = c(2 + 0.01, df_max_mor), x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 0)
  opt_df_zm_global_mor <- optimise(f = zeta_loss, interval = c(2, df_max_mor_global), x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 1)
  opt_df_zd_global_mor <- optimise(f = zeta_loss_deriv, interval = c(2 + 0.01, df_max_mor_global), x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 0)

  # Extract optimal degrees of freedom based on constrained and unconstrained zeta loss and zeta deriv loss using AR data
  opt_df_zm_ar <- optimise(f = zeta_loss, interval = c(2, df_max_ar), x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 1) #WARNS
  opt_df_zd_ar <- optimise(f = zeta_loss_deriv, interval = c(2 + 0.01, df_max_ar), x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 0) #WARNS
  opt_df_zm_global_ar <- optimise(f = zeta_loss, interval = c(2, df_max_ar_global), x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 1) #WARNS
  opt_df_zd_global_ar <- optimise(f = zeta_loss_deriv, interval = c(2 + 0.01, df_max_ar_global), x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 0) #WARNS

  # Check if global search resulted in a considerably better choice of df
  if(opt_df_zm_mor$objective - opt_df_zm_global_mor$objective > 1){
    opt_df_zm_mor <- opt_df_zm_global_mor
  }
  if(opt_df_zd_mor$objective - opt_df_zd_global_mor$objective > 0.5){
    opt_df_zd_mor <- opt_df_zd_global_mor
  }
  if(opt_df_zm_ar$objective - opt_df_zm_global_ar$objective > 1){
    opt_df_zm_ar <- opt_df_zm_global_ar
  }
  if(opt_df_zd_ar$objective - opt_df_zd_global_ar$objective > 0.5){
    opt_df_zd_ar <- opt_df_zd_global_ar
  }

  # Calculate non-linearity scores
  area_df2_df3_mor <- tryCatch(expr = integrate(zeta_loss, lower = 2, upper = 3, x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 0),
                               error = function(e) list(value = NA))
  area_df2_df3_ar <- tryCatch(expr = integrate(zeta_loss, lower = 2, upper = 3, x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 0),
                              error = function(e) list(value = NA))
  darea_df2_df3_mor <- list(value = NA)
  darea_df2_df3_ar <- list(value = NA)

  candidate_df <- seq(from = 2, to = 3, length.out = 2e2L)
  curve_df2_df3_mor <- zeta_loss(df = candidate_df, x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, target = 0)
  curve_df2_df3_ar <- zeta_loss(df = candidate_df, x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, target = 0)
  darea_df2_df3_mor <- list(value = -mean(diff(curve_df2_df3_mor, lag = 2) / 2 / (candidate_df[2] - candidate_df[1])))
  darea_df2_df3_ar <- list(value = -mean(diff(curve_df2_df3_ar, lag = 2) / 2 / (candidate_df[2] - candidate_df[1])))

  # Check if df > 2 is probable to be wrong...
  if((!is.na(area_df2_df3_mor$value) & area_df2_df3_mor$value - 1 <= 1) | (!is.na(darea_df2_df3_mor$value) & darea_df2_df3_mor$value <= 1)){
    if((!is.na(area_df2_df3_mor$value) & area_df2_df3_mor$value - 1 > 1) & (!is.na(darea_df2_df3_mor$value) & darea_df2_df3_mor$value <= 1)){
      opt_df_zm_mor$minimum <- 2
      opt_df_zd_mor$minimum <- 2
    }
    else if((!is.na(area_df2_df3_mor$value) & area_df2_df3_mor$value - 1 <= 1) & (!is.na(darea_df2_df3_mor$value) & darea_df2_df3_mor$value <= 1)){
      opt_df_zm_mor$minimum <- 2
      opt_df_zd_mor$minimum <- 2
    }

  }
  if((!is.na(area_df2_df3_ar$value) & area_df2_df3_ar$value - 1 <= 1 / R) | (!is.na(darea_df2_df3_ar$value) & darea_df2_df3_ar$value <= 1 / R)){
    if((!is.na(area_df2_df3_ar$value) & area_df2_df3_ar$value - 1 > 1 / R) & (!is.na(darea_df2_df3_ar$value) & darea_df2_df3_ar$value <= 1 / R)){
      opt_df_zm_ar$minimum <- 2
      opt_df_zd_ar$minimum <- 2
    }
    else if((!is.na(area_df2_df3_ar$value) & area_df2_df3_ar$value - 1 <= 1 / R) & (!is.na(darea_df2_df3_ar$value) & darea_df2_df3_ar$value <= 1 / R)){
      opt_df_zm_ar$minimum <- 2
      opt_df_zd_ar$minimum <- 2
    }

  }

  #cat("\nNon-linearity score (1MOR):", area_df2_df3_mor$value - 1,
  #    "\nNon-linearity score (2MOR):", darea_df2_df3_mor$value,
  #    "\nNon-linearity score (1AR):", area_df2_df3_ar$value - 1,
  #    "\nNon-linearity score (2AR):", darea_df2_df3_ar$value)

  # Fit smoothing spline objects with respective optimal values of df
  ss_loo <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, df = ss_loo$df, all.knots = all_knots_mor, cv = FALSE)
  ss_gcv <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, df = ss_gcv$df, all.knots = all_knots_mor, cv = FALSE)
  ss_zm_mor <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, df = opt_df_zm_mor$minimum, all.knots = all_knots_mor, cv = FALSE)
  ss_zd_mor <- smooth.spline(x = x_mor, y = y_mor, w = weight_data_mor, df = opt_df_zd_mor$minimum, all.knots = all_knots_mor, cv = FALSE)
  ss_zm_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = opt_df_zm_ar$minimum, all.knots = all_knots_ar, cv = FALSE)
  ss_zd_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = opt_df_zd_ar$minimum, all.knots = all_knots_ar, cv = FALSE)

  ss_loo_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = ss_loo$df, all.knots = all_knots_ar, cv = FALSE)
  ss_gcv_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = ss_gcv$df, all.knots = all_knots_ar, cv = FALSE)
  ss_zm_mor_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = opt_df_zm_mor$minimum, all.knots = all_knots_ar, cv = FALSE)
  ss_zd_mor_ar <- smooth.spline(x = x, y = y, w = weight_data_ar, df = opt_df_zd_mor$minimum, all.knots = all_knots_ar, cv = FALSE)

  # Gather into a list
  ss_list <- list("zm_ar" = ss_zm_ar,
                  "zd_ar" = ss_zd_ar,
                  "zm_mor" = ss_zm_mor_ar,
                  "zd_mor" = ss_zd_mor_ar,
                  "loo" = ss_loo_ar,
                  "gcv" = ss_gcv_ar)

  zeta_list <- lapply(X = ss_list, FUN = function(ss_obj){
    ss_opt_mean_squared_slope <- mean(predict(ss_obj, x = ss_obj$x, deriv = 1)$y^2, na.rm = TRUE)
    var_eps_opt <- ss_obj$pen.crit / (ss_obj$n - ss_obj$df)
    multiplier <- ifelse(ss_obj$n / ss_zm_ar$n < 1/(R-1), 1/R, 1)
    zeta_opt <- var_eps_opt / (var_v * multiplier + var_h * multiplier * ss_opt_mean_squared_slope)
    return(round(zeta_opt, 3L))
  })

  ss_list <- list("zm_ar" = ss_zm_ar,
                  "zd_ar" = ss_zd_ar,
                  "zm_mor" = ss_zm_mor,
                  "zd_mor" = ss_zd_mor,
                  "loo" = ss_loo,
                  "gcv" = ss_gcv)

  zeta_list <- c(list("par" = "zeta"), zeta_list)

  df_list <- list("zm_ar" = opt_df_zm_ar$minimum,
                  "zd_ar" = opt_df_zd_ar$minimum,
                  "zm_mor" = opt_df_zm_mor$minimum,
                  "zd_mor" = opt_df_zd_mor$minimum,
                  "loo" = ifelse(zeta_list$loo <= 0.776, NA, opt_df_loo),
                  "gcv" = ifelse(zeta_list$gcv <= 0.776, NA, opt_df_gcv))
  df_list <- lapply(df_list, round, digits = 3L)

  df_list <- c(list("par" = "df"), df_list)

  residual_scatter_plots <- NULL
  residual_histogram_plots <- NULL
  zeta_vs_df_plots_mor <- NULL
  zeta_vs_df_plots_ar <- NULL

  output_table <- NULL
  output_plots <- NULL

  if(ns_score){
    score_list <- list("zm_ar" = area_df2_df3_ar$value - 1 ,
                       "zd_ar" = darea_df2_df3_ar$value,
                       "zm_mor" = area_df2_df3_mor$value - 1 ,
                       "zd_mor" = darea_df2_df3_mor$value,
                       "loo" = NA,
                       "gcv" = NA)
    score_list <- lapply(score_list, round, digits = 3L)
    score_list <- c(list("par" = "non-linearity score"), score_list)
    output_table <- rbindlist(list(df_list, zeta_list, score_list))
  }
  else{
    output_table <- rbindlist(list(df_list, zeta_list))
  }

  if(plots){
    residual_scatter_plots <- lapply(X = ss_list, FUN = smoothing_spline_residual_plots, which_plot = "scatter_studentized")
    residual_histogram_plots <- lapply(X = ss_list, FUN = smoothing_spline_residual_plots, which_plot = "histogram_studentized")
    zeta_vs_df_plots_mor <- zeta_df_plots(x = x_mor, y = y_mor, weights = weight_data_mor, var_v = var_v, var_h = var_h, all_knots = all_knots_mor, loss = "abs", R = R, mor = TRUE, which_plot = "all")
    zeta_vs_df_plots_ar <- zeta_df_plots(x = x, y = y, weights = weight_data_ar, var_v = var_v, var_h = var_h, all_knots = all_knots_ar, loss = "abs", R = R, mor = FALSE, which_plot = "all")
    output_plots <- list("residual_scatter_plots" = residual_scatter_plots,
                         "residual_histogram_plots" = residual_histogram_plots,
                         "zeta_vs_df_plots_mor" = zeta_vs_df_plots_mor,
                         "zeta_vs_df_plots_ar" = zeta_vs_df_plots_ar,
                         "second_deriv" = second_deriv_object)
    output <- list("table" = output_table, "plots" = output_plots)
    if(models){
      output <- c(output, list("models" = ss_list))
    }

  }
  else{
    output <- output_table
    if(models){
      output <- c(list("table" = output), list("models" = ss_list))
    }
  }

  return(output)
}

#' Comparison-wise Suggestions for Optimal Starting Values for Effective Degrees of Freedom
#'
#' @param data A \code{list} or \code{data.table} containing columns \code{comparison},
#'             \code{SampleID}, \code{ReplicateID}, \code{MP_A} and \code{MP_B}.
#' @param ns_score A \code{logical}. If \code{TRUE}, calculates non-linearity scores for
#'                 each method used to derive optimal effective DOF values. See details.
#' @param use_weights A code{logical}. If \code{TRUE}, uses weights in the algorithms.
#'                    Weights are automatically determined based on \code{data}.
#' @param output_type A \code{character} value.
#' @param na_rm A \code{logical}. If \code{TRUE}, removes NA values before modeling. Note
#'              that this affects algorithm results but not measurement variance or replicate
#'              count calculations.
#'
#' @return A \code{data.table}
#' @export
#'
#' @examples print(1)
suggest_dfs <- function(data, ns_score = FALSE, use_weights = FALSE, output_type = c("mean", "median", "zd", "zm"), na_rm = TRUE){
  suggest_df_data <- data[, suggest_df(.SD, ns_score = ns_score, use_weights = use_weights, na_rm = na_rm), by = comparison][par == "df"]
  output_type <- output_type[1]
  if(output_type == "mean"){
    suggest_df_data <- melt.data.table(suggest_df_data, id.col = "comparison", measure.vars = 3:7, value.name = "df")[, list(df = mean(df, na.rm = TRUE)), by = comparison]
  }
  else if(output_type == "median"){
    suggest_df_data <- melt.data.table(suggest_df_data, id.col = "comparison", measure.vars = 3:7, value.name = "df")[, list(df = median(df, na.rm = TRUE)), by = comparison]
  }
  else if(output_type == "zd"){
    suggest_df_data <- suggest_df_data[, c("comparison", "zd_mor")]
    names(suggest_df_data)[which("zd_mor" == names(suggest_df_data))] <- "df"
  }
  else if(output_type == "zm"){
    suggest_df_data <- suggest_df_data[, c("comparison", "zm_mor")]
    names(suggest_df_data)[which("zm_mor" == names(suggest_df_data))] <- "df"
  }

  estimated_zetas <- estimate_zetas_ss(data = data, df = suggest_df_data$df, use_weights = use_weights, na_rm = na_rm)
  return(merge(suggest_df_data, estimated_zetas, by = "comparison"))
}

