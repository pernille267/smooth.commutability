#' Internal Function for Constructing List of Additional Parameters
#'
#' This function is used internally to readily get additional_parameters.
#' @keywords internal
get_additional_parameters <- function(smoothing_type = "loocv", smudge = 1.0){
  if(smoothing_type == "loocv"){
    return(list(method = "loocv", smudge = smudge, c_star_min = 0.00, c_star_max = 1.50, tol = 0.5, window = 10, iter = 1))
  }
  else if(smoothing_type == "loocv_restricted"){
    return(list(method = "loocv", smudge = smudge, c_star_min = 0.38, c_star_max = 1.50, tol = 0.5, window = 10, iter = 1))
  }
  else if(smoothing_type == "gcv"){
    return(list(method = "gcv", smudge = smudge, c_star_min = 0.00, c_star_max = 1.50, tol = 0.5, window = 10, iter = 1))
  }
  else if(smoothing_type == "gcv_restricted"){
    return(list(method = "gcv", smudge = smudge, c_star_min = 0.38, c_star_max = 1.50, tol = 0.5, window = 10, iter = 1))
  }
  else{
    stop("ERROR")
  }
}

#' Internal Function for Estimating Optimal Values of DOF
#'
#' This function is used internally to readily get optimal values of \eqn{\mathrm{df}(c)}.
#' @keywords internal
optimal_df <- function(data){
  df_loocv <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv"))
  df_loocv_restricted <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv_restricted"))
  df_loocv_weighted <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv"))
  df_loocv_restricted_weighted <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv_restricted"))

  df_gcv <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv"))
  df_gcv_restricted <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted"))
  df_gcv_weighted <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv"))
  df_gcv_restricted_weighted <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted"))

  df_gcv_smudge <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv", 1.2))
  df_gcv_restricted_smudge <- smoothing_spline(data = data, weights = 1, df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted", 1.2))
  df_gcv_weighted_smudge <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv", 1.2))
  df_gcv_restricted_weighted_smudge <- smoothing_spline(data = data, weights = "estimate", df_max = NULL, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted", 1.2))

  optimal_dfs <- c(df_loocv$df,
                   df_loocv_restricted$df,
                   df_loocv_weighted$df,
                   df_loocv_restricted_weighted$df,
                   df_gcv$df,
                   df_gcv_restricted$df,
                   df_gcv_weighted$df,
                   df_gcv_restricted_weighted$df,
                   df_gcv_smudge$df,
                   df_gcv_restricted_smudge$df,
                   df_gcv_weighted_smudge$df,
                   df_gcv_restricted_weighted_smudge$df)

  cv_crits <- c(df_loocv$cv_crit,
                df_loocv_restricted$cv_crit,
                df_loocv_weighted$cv_crit,
                df_loocv_restricted_weighted$cv_crit,
                df_gcv$cv_crit,
                df_gcv_restricted$cv_crit,
                df_gcv_weighted$cv_crit,
                df_gcv_restricted_weighted$cv_crit,
                df_gcv_smudge$cv_crit,
                df_gcv_restricted_smudge$cv_crit,
                df_gcv_weighted_smudge$cv_crit,
                df_gcv_restricted_weighted_smudge$cv_crit)

  general_method <- c(rep("LOOCV", 4L), rep("GCV", 8L))
  restricted <- rep(c("no", "yes", "no", "yes"), 3L)
  weighted <- rep(c("no", "no", "yes", "yes"), 3L)
  smudged <- c(rep("no", 8L), rep("yes", 4L))

  optimal_dfs <- data.table("df" = optimal_dfs,
                            "method" = general_method,
                            "restricted" = restricted,
                            "weighted" = weighted,
                            "smudged" = smudged,
                            "cv" = cv_crits)

  return(optimal_dfs)

}

#' Internal Function for Estimating LOO and Generalized Cross-Validation Scores
#'
#' This function is used internally to readily get cross-validation scores.
#' @keywords internal
calculate_cv_score <- function(data, df = 2){
  cv_crit_loocv <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv"))$cv_crit
  cv_crit_loocv_restricted <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv_restricted"))$cv_crit
  cv_crit_loocv_weighted <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv"))$cv_crit
  cv_crit_loocv_restricted_weighted <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("loocv_restricted"))$cv_crit

  cv_crit_gcv <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv"))$cv_crit
  cv_crit_gcv_restricted <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted"))$cv_crit
  cv_crit_gcv_weighted <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv"))$cv_crit
  cv_crit_gcv_restricted_weighted <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted"))$cv_crit

  cv_crit_gcv_smudge <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv", 1.2))$cv_crit
  cv_crit_gcv_restricted_smudge <- smoothing_spline(data = data, weights = 1, df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted", 1.2))$cv_crit
  cv_crit_gcv_weighted_smudge <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv", 1.2))$cv_crit
  cv_crit_gcv_restricted_weighted_smudge <- smoothing_spline(data = data, weights = "estimate", df = df, attempt_fast = FALSE,  additional_parameters = get_additional_parameters("gcv_restricted", 1.2))$cv_crit

  cv_score <- c(cv_crit_loocv,
                cv_crit_loocv_restricted,
                cv_crit_loocv_weighted,
                cv_crit_loocv_restricted_weighted,
                cv_crit_gcv,
                cv_crit_gcv_restricted,
                cv_crit_gcv_weighted,
                cv_crit_gcv_restricted_weighted,
                cv_crit_gcv_smudge,
                cv_crit_gcv_restricted_smudge,
                cv_crit_gcv_weighted_smudge,
                cv_crit_gcv_restricted_weighted_smudge)

  general_method <- c(rep("LOOCV", 4L), rep("GCV", 8L))
  restricted <- rep(c("no", "yes", "no", "yes"), 3L)
  weighted <- rep(c("no", "no", "yes", "yes"), 3L)
  smudged <- c(rep("no", 8L), rep("yes", 4L))

  cv_score <- data.table("df" = rep(df, 12L),
                         "method" = general_method,
                         "restricted" = restricted,
                         "weighted" = weighted,
                         "smudged" = smudged,
                         "cv" = cv_score)

  return(cv_score)
}

#' Internal Function for Calculating \eqn{\mathrm{df}(c)} Versus Cross-Validation Scores
#'
#' This function is used internally to calculate cross-validation scores for a grid of \eqn{\mathrm{df}(c)} values.
#' @keywords internal
df_vs_cv <- function(data, df_grid){
  cv_scores <- sapply(X = df_grid, FUN = function(x){
    calculate_cv_score(data, x)
  }, simplify = FALSE)
  cv_scores <- rbindlist(cv_scores)
  return(cv_scores)
}

#' Internal Function for Rendering \eqn{\mathrm{df}(c)} Versus Cross-Validation Score Plots
#'
#' This function is used internally to plot \eqn{\mathrm{df}(c)} against cross-validation scores.
#' @keywords internal
df_vs_cv_plot <- function(plot_data, n_breaks_y = 10, log_scale = FALSE, method_colors = c("black", "red"), smudge_linetypes = c("longdash", "solid")){

  cv <- method <- smudged <- weighted <- restricted <- NULL

  # Relevant plotting components
  df_breaks <- seq(from = 2,
                   to = max(plot_data$df, na.rm = TRUE),
                   by = 1)

  trans_y <- if(log_scale){"log10"}else{"identity"}
  color_scale <- c("LOOCV" = method_colors[1],
                   "GCV" = method_colors[2])
  linetype_scale <- c("yes" = smudge_linetypes[1],
                      "no" = smudge_linetypes[2])

  plot_output <- ggplot(data = plot_data) +
    geom_line(mapping = aes(x = df, y = cv, color = method, linetype = smudged)) +
    scale_x_continuous(name = "df(c)",
                       breaks = df_breaks,
                       expand = c(0.02, 0.02)) +
    scale_y_continuous(name = "Cross-Validation Score",
                       n.breaks = n_breaks_y,
                       transform = trans_y,
                       labels = function(x) round(x, 3L)) +
    scale_color_manual(values = color_scale) +
    scale_linetype_manual(values = linetype_scale,
                          labels = c("1", "1.2")) +
    labs(color = "Method:",
         linetype = "Smudge Factor:")

  if((!is.null(plot_data$weighted)) & (!is.null(plot_data$restricted))){
    plot_output <- plot_output +
      facet_grid(rows = vars(weighted),
                 cols = vars(restricted),
                 labeller = labeller(.rows = function(x) paste0("Weighted: ", x),
                                     .cols = function(x) paste0("Restricted: ", x)),
                 scales = "free_y") +
      theme_bw() +
      theme(axis.title = element_text(color = "#000000"),
            axis.text = element_text(color = "#000000"),
            legend.background = element_rect(fill = "#FFFFFF", color = "black"),
            legend.key = element_rect(fill = "#FFFFFF", color = "black"),
            legend.title = element_text(face = "bold"),
            legend.position = "top",
            strip.background = element_rect(fill = "#000000"),
            strip.text = element_text(face = "bold", color = "#FFFFFF"))
  }
  else if(is.null(plot_data$weighted) & (!is.null(plot_data$restricted))){
    plot_output <- plot_output +
      facet_grid(cols = vars(restricted),
                 labeller = labeller(.cols = function(x) paste0("Restricted: ", x))) +
      theme_bw() +
      theme(axis.title = element_text(color = "#000000"),
            axis.text = element_text(color = "#000000"),
            legend.background = element_rect(fill = "#FFFFFF", color = "black"),
            legend.key = element_rect(fill = "#FFFFFF", color = "black"),
            legend.title = element_text(face = "bold"),
            legend.position = "top",
            strip.background = element_rect(fill = "#000000"),
            strip.text = element_text(face = "bold", color = "#FFFFFF"))
  }
  else if((!is.null(plot_data$weighted)) & is.null(plot_data$restricted)){
    plot_output <- plot_output +
      facet_grid(cols = vars(weighted),
                 labeller = labeller(.cols = function(x) paste0("Weighted: ", x)),
                 scales = "free") +
      theme_bw() +
      theme(axis.title = element_text(color = "#000000"),
            axis.text = element_text(color = "#000000"),
            legend.background = element_rect(fill = "#FFFFFF", color = "black"),
            legend.key = element_rect(fill = "#FFFFFF", color = "black"),
            legend.title = element_text(face = "bold"),
            legend.position = "top",
            strip.background = element_rect(fill = "#000000"),
            strip.text = element_text(face = "bold", color = "#FFFFFF"))
  }
  else{
    plot_output <- plot_output +
      theme_bw() +
      theme(axis.title = element_text(color = "#000000"),
            axis.text = element_text(color = "#000000"),
            legend.background = element_rect(fill = "#FFFFFF", color = "black"),
            legend.key = element_rect(fill = "#FFFFFF", color = "black"),
            legend.title = element_text(face = "bold"),
            legend.position = "top",
            strip.background = element_rect(fill = "#000000"),
            strip.text = element_text(face = "bold", color = "#FFFFFF"))
  }

  return(plot_output)
}

#' Internal Function for Calculating \eqn{\hat{\zeta}} for each \eqn{\mathrm{df}(c)}
#'
#' This function is used internally to calculate \eqn{\hat{\zeta}} for each \eqn{\mathrm{df}(c)}.
#' @keywords internal
df_vs_zeta <- function(data, df_grid){

  # Calculate zeta hat values for each df along the grid.
  zeta_hats_ar <- sapply(X = df_grid,
                         FUN = function(x){
                           estimate_zeta_ss(data = data,
                                            df = x,
                                            weighted = FALSE,
                                            mor = FALSE,
                                            na_rm = TRUE)$zeta
                         })
  zeta_hats_ar_weighted <- sapply(X = df_grid,
                                  FUN = function(x){
                                    estimate_zeta_ss(data = data,
                                                     df = x,
                                                     weighted = TRUE,
                                                     mor = FALSE,
                                                     na_rm = TRUE)$zeta
                                  })
  zeta_hats_mor <- sapply(X = df_grid,
                          FUN = function(x){
                            estimate_zeta_ss(data = data,
                                             df = x,
                                             weighted = FALSE,
                                             mor = TRUE,
                                             na_rm = TRUE)$zeta
                          })
  zeta_hats_mor_weighted <- sapply(X = df_grid,
                                   FUN = function(x){
                                     estimate_zeta_ss(data = data,
                                                      df = x,
                                                      weighted = TRUE,
                                                      mor = TRUE,
                                                      na_rm = TRUE)$zeta
                                   })

  # Combine results into a data.table..
  zeta_hats <- data.table(df = rep(df_grid, 4L),
                          mor = rep(c("no", "yes"),
                                    each = 2L * length(df_grid)),
                          weighted = rep(rep(c("no", "yes"),
                                             each = length(df_grid)),
                                         2L),
                          zeta = c(zeta_hats_ar,
                                   zeta_hats_ar_weighted,
                                   zeta_hats_mor,
                                   zeta_hats_mor_weighted))
  return(zeta_hats)
}

#' Internal Function for Rendering \eqn{\mathrm{df}(c)} Versus \eqn{\hat{\zeta}} Plots
#'
#' This function is used internally to plot \eqn{\mathrm{df}(c)} against \eqn{\hat{\zeta}}.
#' @keywords internal
df_vs_zeta_plot <- function(plot_data, n_breaks_y = 10, log_scale = FALSE, data_structure_colors = c("black", "red"), weighted_linetypes = c("longdash", "solid")){

  zeta <- mor <- weighted <- NULL

  # Relevant plotting components
  df_breaks <- seq(from = 2,
                   to = max(plot_data$df, na.rm = TRUE),
                   by = 1)

  trans_y <- if(log_scale){"log10"}else{"identity"}
  color_scale <- c("yes" = data_structure_colors[1],
                   "no" = data_structure_colors[2])
  linetype_scale <- c("yes" = weighted_linetypes[1],
                      "no" = weighted_linetypes[2])


  plot_output <- ggplot(data = plot_data) +
    scale_x_continuous(name = "df(c)",
                       breaks = df_breaks,
                       expand = c(0.02, 0.02)) +
    scale_y_continuous(name = expression(hat(zeta)),
                       n.breaks = n_breaks_y,
                       transform = trans_y,
                       labels = function(x) round(x, 3L))


  if(!is.null(plot_data$weighted)){
    plot_output <- plot_output +
      geom_line(mapping = aes(x = df, y = zeta, color = mor, linetype = weighted)) +
      scale_color_manual(values = color_scale,
                         labels = c("AR", "MOR")) +
      scale_linetype_manual(values = linetype_scale,
                            labels = c("No", "Yes")) +
      labs(color = "Data Structure: ",
           linetype = "Weighted: ")
  }
  else{
    plot_output <- plot_output +
      geom_line(mapping = aes(x = df, y = zeta, color = mor)) +
      scale_color_manual(values = color_scale,
                         labels = c("AR", "MOR")) +
      labs(color = "Data Structure: ")
  }

  plot_output <- plot_output +
    theme_bw() +
    theme(axis.title.y = element_text(color = "#000000", angle = 0, vjust = 0.5),
          axis.title.x = element_text(color = "#000000"),
          axis.text = element_text(color = "#000000"),
          legend.background = element_rect(fill = "#FFFFFF", color = "black"),
          legend.key = element_rect(fill = "#FFFFFF", color = "black"),
          legend.title = element_text(face = "bold"),
          legend.position = "top")

  return(plot_output)

}


#' Smoothing Spline Diagnostics
#'
#' @param data A dataset
#' @param only_optimal_dfs A \code{logical} value. If \code{TRUE}, only the
#'                         optimal values of \eqn{\mathrm{df}(c)} are returned.
#'                         Not plots and no tables contaning relationships between
#'                         \eqn{\mathrm{df}(c)} and \eqn{\hat{\zeta}} and cross-validation
#'                         scores.
#' @param weighted A \code{logical} value. If \code{TRUE}, include weighted results.
#' @param exclude_not_restricted A \code{logical} value. If \code{FALSE}, exclude
#'                               non-restricted results.
#' @param na_rm A \code{logical} value. If \code{TRUE}, \code{NA} values are
#'        removed before the diagnostics are performed.
#'
#' @description
#' Perform diagnostics of a smoothing spline fit based.
#'
#' @details
#' Before utilizing the \code{smoothing_spline()} function to fit a smoothing spline
#' to the \code{data}, it may be advantegous to consider some diagnostics of the
#' regarded dataset. This function return key diagnostic components such as
#' \itemize{
#'  \item \eqn{\mathrm{df}(c)} versus Cross-Validation scores
#'  \item \eqn{\mathrm{df}(c)} versus \eqn{\hat{\zeta}}
#'  \item Optimal values of \eqn{\mathrm{df}(c)} based on different variants of
#'        Generalized Cross-Validation (GCV) and Leave-One-Out (LOO) Cross-Validation
#' }
#' These diagnostic components may aid in choosing favorable values of \eqn{\mathrm{df}(c)}
#' based on \code{data}.
#'
#' @return
#' Returns a \code{list}:
#' \itemize{
#'  \item \code{item1:} bla bla
#'  \item \code{item2:} bla bla
#'  \item \code{item3:} bla bla
#' }
#'
#' @export
#'
#' @examples
#' print(1)
smoothing_spline_diagnostics <- function(data, only_optimal_dfs = FALSE, weighted = FALSE, exclude_not_restricted = TRUE, na_rm = TRUE){

  restricted <- NULL
  df_grid <- seq(from = 2, to = 10, by = 0.25)

  # Check for required variables in data
  if(!any("MP_B" == names(data)) | !any("MP_A" == names(data)) | !any("SampleID" == names(data)) | !any("ReplicateID" == names(data))){
    stop("'data' is expected to have variables named 'MP_B' and 'MP_A'.")
  }

  # Calculate imprecision estimates
  impr <- global_precision_estimates(data)

  # Calculate MOR data
  mor_data <- fun_of_replicates2(data)

  # Remove NA-values if na_rm = TRUE
  if(isTRUE(na_rm)){
    if(any(is.na(data$MP_B)) | any(is.na(data$MP_A))){
      data <- na.omit(data)
    }
  }

  if(only_optimal_dfs){

    optimal_dfs <- optimal_df(mor_data)

    # Remove weighted
    if(!weighted){
      optimal_dfs <- optimal_dfs[weighted == "no"]
      optimal_dfs$weighted <- NULL
    }

    # Remove non-restricted
    if(exclude_not_restricted){
      optimal_dfs <- optimal_dfs[restricted == "yes"]
      optimal_dfs$restricted <- NULL
    }

    return(optimal_dfs)

  }

  # Get components
  df_zeta <- df_vs_zeta(data, df_grid)
  df_cv <- df_vs_cv(mor_data, df_grid)


  # Remove weighted
  if(!weighted){
    df_zeta <- df_zeta[weighted == "no"]
    df_cv <- df_cv[weighted == "no"]
    optimal_dfs <- optimal_dfs[weighted == "no"]
    df_zeta$weighted <- NULL
    df_cv$weighted <- NULL
    optimal_dfs$weighted <- NULL
  }

  # Remove non-restricted
  if(exclude_not_restricted){
    df_cv <- df_cv[restricted == "yes"]
    optimal_dfs <- optimal_dfs[restricted == "yes"]
    df_cv$restricted <- NULL
    optimal_dfs$restricted <- NULL
  }

  # Gather components
  output <- list("optimal_dfs" = optimal_dfs,
                 "df_vs_cv" = df_cv,
                 "df_vs_zeta" = df_zeta,
                 "df_vs_cv_plot" = df_vs_cv_plot(df_cv),
                 "df_vs_zeta_plot" = df_vs_zeta_plot(df_zeta))

  return(output)

}
