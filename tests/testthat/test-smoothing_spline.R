library(fasteqa)
library(data.table)
library(testthat)

# Reproducibility
set.seed(99)

# Test data
test_data_1 <- sim_eqa_data(parameters = list(n = 25,
                                              R = 3,
                                              cvx = 0.01,
                                              cvy = 0.01,
                                              cil = 2,
                                              ciu = 10),
                            type = 2,
                            AR = FALSE) |> setDT()
test_data_2_ar <- fix_predictor_and_response(crp_cs_data)
test_data_2_mor <- test_data_2_ar[, fun_of_replicates(.SD), by = comparison]

# Output from smoothing_spline()
ss_fit_1 <- smoothing_spline(data = test_data_1, df = 3)
ss_fit_2 <- smoothing_spline(data = test_data_1, df = 5)
ss_fit_3 <- smoothing_spline(data = test_data_1, df = 7.5)

# Output from smooth.spline()
ss_fit_1_ref <- smooth.spline(x = test_data_1$MP_B, test_data_1$MP_A, df = 3, keep.stuff = TRUE, all.knots = c(0, ss_fit_1$interior_knots, 1), cv = TRUE)
ss_fit_2_ref <- smooth.spline(x = test_data_1$MP_B, test_data_1$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_2$interior_knots, 1), cv = TRUE)
ss_fit_3_ref <- smooth.spline(x = test_data_1$MP_B, test_data_1$MP_A, df = 7.5, keep.stuff = TRUE, , all.knots = c(0, ss_fit_3$interior_knots, 1), cv = TRUE)

# Reconstruct matrices from smooth.spline()
ss_fit_1_ref_matrices <- get_matrices(ss_fit_1_ref$auxM)
ss_fit_2_ref_matrices <- get_matrices(ss_fit_2_ref$auxM)
ss_fit_3_ref_matrices <- get_matrices(ss_fit_3_ref$auxM)


test_that("Check if df corresponds with the given", {
  expect_true(object = abs(ss_fit_1$df - 3) < 1e-3)
  expect_true(object = abs(ss_fit_2$df - 5) < 1e-3)
  expect_true(object = abs(ss_fit_3$df - 7.5) < 1e-3)
})

test_that("Check if lambda corresponds with those given by smooth.spline", {
  expect_true(object = abs(ss_fit_1$lambda - ss_fit_1_ref$lambda) < 1e-3)
  expect_true(object = abs(ss_fit_2$lambda - ss_fit_2_ref$lambda) < 1e-3)
  expect_true(object = abs(ss_fit_3$lambda - ss_fit_3_ref$lambda) < 1e-3)
})

test_that("Check if variance estimators corresponds with those given by smooth.spline", {
  expect_true(object = abs(ss_fit_1$var_eps - ss_fit_1_ref$pen.crit / (25 - ss_fit_1_ref$df)) < 1e-3)
  expect_true(object = abs(ss_fit_2$var_eps - ss_fit_2_ref$pen.crit / (25 - ss_fit_2_ref$df)) < 1e-3)
  expect_true(object = abs(ss_fit_3$var_eps - ss_fit_3_ref$pen.crit / (25 - ss_fit_3_ref$df)) < 1e-3)
})

test_that("Check if fitted values correspond with those given by smooth.spline", {
  expect_true(object = all(abs(ss_fit_1$fitted[,] - ss_fit_1_ref$y) < 1e-3))
  expect_true(object = all(abs(ss_fit_2$fitted[,] - ss_fit_2_ref$y) < 1e-3))
  expect_true(object = all(abs(ss_fit_3$fitted[,] - ss_fit_3_ref$y) < 1e-3))
})

test_that("Check if leverage values correspond with those given by smooth.spline", {
  expect_true(object = all(abs(diag(ss_fit_1$smoothing_matrix) - ss_fit_1_ref$lev) < 1e-3))
  expect_true(object = all(abs(diag(ss_fit_2$smoothing_matrix) - ss_fit_2_ref$lev) < 1e-3))
  expect_true(object = all(abs(diag(ss_fit_3$smoothing_matrix) - ss_fit_3_ref$lev) < 1e-3))
})

test_that("Check if penalty matrix is approximately equal with that by smooth.spline", {
  expect_true(median(abs((ss_fit_1$penalty_matrix - ss_fit_1_ref_matrices$Omega) / ss_fit_1_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
  expect_true(median(abs((ss_fit_2$penalty_matrix - ss_fit_2_ref_matrices$Omega) / ss_fit_2_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
  expect_true(median(abs((ss_fit_3$penalty_matrix - ss_fit_3_ref_matrices$Omega) / ss_fit_3_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
})

test_that("Check if penalty matrix is approximately equal with that by smooth.spline", {
  expect_true(median(abs((ss_fit_1$penalty_matrix - ss_fit_1_ref_matrices$Omega) / ss_fit_1_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
  expect_true(median(abs((ss_fit_2$penalty_matrix - ss_fit_2_ref_matrices$Omega) / ss_fit_2_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
  expect_true(median(abs((ss_fit_3$penalty_matrix - ss_fit_3_ref_matrices$Omega) / ss_fit_3_ref_matrices$Omega) * 100, na.rm = TRUE) < 0.5)
})

Q1 <- solve(ss_fit_1_ref_matrices$BTB + ss_fit_1_ref_matrices$Omega * ss_fit_1_ref$lambda)
Q2 <- solve(ss_fit_2_ref_matrices$BTB + ss_fit_2_ref_matrices$Omega * ss_fit_2_ref$lambda)
Q3 <- solve(ss_fit_3_ref_matrices$BTB + ss_fit_3_ref_matrices$Omega * ss_fit_3_ref$lambda)

var_beta_ref_1 <- ss_fit_1_ref$pen.crit / (25 - ss_fit_1_ref$df) * Q1 %*% ss_fit_1_ref_matrices$BTB %*% Q1
var_beta_ref_2 <- ss_fit_2_ref$pen.crit / (25 - ss_fit_2_ref$df) * Q2 %*% ss_fit_2_ref_matrices$BTB %*% Q2
var_beta_ref_3 <- ss_fit_3_ref$pen.crit / (25 - ss_fit_3_ref$df) * Q3 %*% ss_fit_3_ref_matrices$BTB %*% Q3

test_that("Check if covariance matrix of coefficients is approximately equal with that by smooth.spline", {
  expect_true(median(abs((var_beta_ref_1 / ss_fit_1$cov_beta - 1) * 100), na.rm = TRUE) < 0.5)
  expect_true(median(abs((var_beta_ref_2 / ss_fit_2$cov_beta - 1) * 100), na.rm = TRUE) < 0.5)
  expect_true(median(abs((var_beta_ref_3 / ss_fit_3$cov_beta - 1) * 100), na.rm = TRUE) < 0.5)
})

test_that("Check if coefficients correspond with those given by smooth.spline", {
  expect_true(object = all(abs(ss_fit_1$coefficients[,] - ss_fit_1_ref$fit$coef) < 1e-3))
  expect_true(object = all(abs(ss_fit_2$coefficients[,] - ss_fit_2_ref$fit$coef) < 1e-3))
  expect_true(object = all(abs(ss_fit_3$coefficients[,] - ss_fit_3_ref$fit$coef) < 1e-3))
})

test_that("Check if predictions align with ols predictions if df = 2", {
  test_data <- test_data_2_mor[comparison == "AQT90 - QuickRead"]
  impr_data <- global_precision_estimates(test_data_2_ar[comparison == "AQT90 - QuickRead"])
  smooth_spline_fit_1 <- smoothing_spline(data = test_data, df = 2)$fitted[,]
  ols_fit_1 <- predict_eqa(data = test_data, new_data = test_data, imprecision_estimates = impr_data, rounding = 6, method = "ols")$prediction |> sort()
  expect_true(all(abs(smooth_spline_fit_1 - ols_fit_1) < 1e-3))
})

multiple_ss_fits <- lapply(X = split(test_data_2_mor, by = "comparison", keep.by = TRUE), FUN = smoothing_spline)

test_that("Check if printing function works as expected", {
  choices <- c("Smoothing Spline Fit:", "Scale-invariant Smoothing", "Computational lambda",
               "Maximum Degrees of", "Effective Degrees of F", "Estimated Residual Variance: ",
               "LOOCV:")
  lapply(multiple_ss_fits, function(x){
    sample_choice <- sample(choices, size = 1)
    expect_output(object = print.smoothing_spline(x), regexp = sample_choice)
  })
})

test_that("Check if plotting function works as expected", {
  sample_choices <- sample.int(n = 10, size = 3, replace = FALSE)
  output_plot_1 <- plot.smoothing_spline(multiple_ss_fits[[sample_choices[1]]], type = "confidence")
  output_plot_2 <- plot.smoothing_spline(multiple_ss_fits[[sample_choices[2]]], type = "prediction")
  output_plot_3 <- plot.smoothing_spline(multiple_ss_fits[[sample_choices[3]]], type = "none")
  expect_s3_class(object = output_plot_1, class = "ggplot")
  expect_s3_class(object = output_plot_2, class = "ggplot")
  expect_s3_class(object = output_plot_3, class = "ggplot")
  expect_equal(output_plot_1$labels$title, "Smoothing Spline Fit + 95% Confidence Intervals:")
  expect_equal(output_plot_2$labels$title, "Smoothing Spline Fit + 95% Prediction Intervals:")
  expect_equal(output_plot_3$labels$title, "Smoothing Spline Fit:")
})

output1 <- smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = 5, attempt_fast = TRUE)
output2 <- smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = 5, attempt_fast = FALSE)

test_that("Check if attempt_fast gives equivalent results", {
  expect_equal(output1$fitted, output2$fitted[,], tolerance = 1e-4)
  expect_equal(output1$residuals, output2$residuals[,], tolerance = 1e-4)
  expect_equal(output1$df, output2$df, tolerance = 1e-4)
  expect_equal(output1$lambda, output2$lambda, tolerance = 1e-4)
  expect_equal(output1$sp, output2$sp, tolerance = 1e-1)
  #expect_equal(output1$cv_crit, output2$cv_crit, tolerance = 1e-4)
  expect_equal(output1$var_eps, output2$var_eps, tolerance = 1e-4)
  expect_equal(output1$var_fit, output2$var_fit, tolerance = 1e-4)
  expect_equal(output1$coefficients, output2$coefficients[,], tolerance = 1e-4)
  expect_equal(output1$cov_beta, output2$cov_beta, tolerance = 1e-4)
  expect_equal(output1$smoothing_matrix, output2$smoothing_matrix, tolerance = 1e-4)
  expect_equal(output1$lambda * output1$penalty_matrix, output2$lambda * output2$penalty_matrix, tolerance = 1e-1)
  expect_equal(output1$B, output2$B, tolerance = 1e-9)
  expect_equal(output1$BTB, output2$BTB, tolerance = 1e-4)
  expect_equal(output1$interior_knots, output2$interior_knots, tolerance = 1e-4)
  expect_equal(output1$knot_boundary, output2$knot_boundary, tolerance = 1e-9)
  expect_equal(output1$x, output2$x, tolerance = 1e-9)
  expect_equal(output1$y, output2$y, tolerance = 1e-9)
})


test_that("Check if attempt_fast actually is faster", {
  time1_before <- Sys.time()
  output1 <- replicate(n = 5e2, smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = NULL, lambda = NULL, df_max = 25, attempt_fast = TRUE))
  time1_after <- Sys.time()
  time1 <- as.double(time1_after - time1_before)

  time2_before <- Sys.time()
  output2 <- replicate(n = 5e2, smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = NULL, lambda = NULL, df_max = 25, attempt_fast = FALSE))
  time2_after <- Sys.time()
  time2 <- as.double(time2_after - time2_before)
  expect_lte(time1, time2)

  time1_before <- Sys.time()
  output1 <- replicate(n = 5e2, smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = 6, lambda = NULL, df_max = 25, attempt_fast = TRUE))
  time1_after <- Sys.time()
  time1 <- as.double(time1_after - time1_before)

  time2_before <- Sys.time()
  output2 <- replicate(n = 5e2, smoothing_spline(data = split(test_data_2_mor, by = "comparison", keep.by = TRUE)[[1]], df = 6, lambda = NULL, df_max = 25, attempt_fast = FALSE))
  time2_after <- Sys.time()
  time2 <- as.double(time2_after - time2_before)
  expect_lte(time1, time2)
})
