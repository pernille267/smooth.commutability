library(fasteqa)
library(data.table)
library(testthat)

# Reproducibility
set.seed(99)

# Parameters used
cs_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10)
eq_parameters <- c(list(obs_tau = c(1.5, 7, 11)), cs_parameters)

# Simulated clinical sample data
test_data_cs_1 <- sim_eqa_data(parameters = cs_parameters, type = 1) |> setDT()
test_data_cs_2 <- sim_eqa_data(parameters = cs_parameters, type = 2) |> setDT()
test_data_cs_3 <- sim_eqa_data(parameters = cs_parameters, type = 3) |> setDT()

# Simulated evaluated material data
test_data_eq_1 <- sim_eqa_data(parameters = eq_parameters, type = 1) |> setDT()
test_data_eq_2 <- sim_eqa_data(parameters = eq_parameters, type = 2) |> setDT()
test_data_eq_3 <- sim_eqa_data(parameters = eq_parameters, type = 3) |> setDT()
test_data_eq_1[, SampleID := paste0("EQAM_", SampleID)]
test_data_eq_2[, SampleID := paste0("EQAM_", SampleID)]
test_data_eq_3[, SampleID := paste0("EQAM_", SampleID)]

ss_fit_1 <- smoothing_spline(data = test_data_cs_1, df = 5)
ss_fit_2 <- smoothing_spline(data = test_data_cs_2, df = 5)
ss_fit_3 <- smoothing_spline(data = test_data_cs_3, df = 5)
ss_fit_4 <- smoothing_spline(data = test_data_cs_2, df = 2)

ss_fit_1_ref <- smooth.spline(x = test_data_cs_1$MP_B, test_data_cs_1$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_1$interior_knots, 1))
ss_fit_2_ref <- smooth.spline(x = test_data_cs_2$MP_B, test_data_cs_2$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_2$interior_knots, 1))
ss_fit_3_ref <- smooth.spline(x = test_data_cs_3$MP_B, test_data_cs_3$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_3$interior_knots, 1))
ss_fit_4_ref <- smooth.spline(x = test_data_cs_2$MP_B, test_data_cs_2$MP_A, df = 2, keep.stuff = TRUE, all.knots = c(0, ss_fit_2$interior_knots, 1))

ss_pred_1 <- predict_smoothing_spline(data = test_data_cs_1, new_data = test_data_eq_1, df = 5, level = 0.95, rounding = 6)
ss_pred_2 <- predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_2, df = 5, level = 0.95, rounding = 6)
ss_pred_3 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_3, df = 5, level = 0.95, rounding = 6)
ss_pred_4 <- predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_2, df = 2, level = 0.95, rounding = 6)

ss_pred_1_ref <- predict(ss_fit_1_ref, x = test_data_eq_1$MP_B)
ss_pred_2_ref <- predict(ss_fit_2_ref, x = test_data_eq_2$MP_B)
ss_pred_3_ref <- predict(ss_fit_3_ref, x = test_data_eq_3$MP_B)

test_that("Check if predictions correspond with those given by smooth.spline", {
  expect_true(object = all(abs(ss_pred_1$prediction - ss_pred_1_ref$y) < 0.01))
  expect_true(object = all(abs(ss_pred_2$prediction - ss_pred_2_ref$y) < 0.01))
  expect_true(object = all(abs(ss_pred_3$prediction - ss_pred_3_ref$y) < 0.01))
})

test_that("Check if extrapolation is correctly registered", {
  expect_true(object = all(ss_pred_1$extrapolate == c(1, 0, 1)))
  expect_true(object = all(ss_pred_2$extrapolate == c(1, 0, 1)))
  expect_true(object = all(ss_pred_3$extrapolate == c(1, 0, 1)))
})

ols_pred_1 <- predict_eqa(data = test_data_cs_1, new_data = test_data_eq_1, method = "ols", level = 0.95, imprecision_estimates = list(Var_A = 3, Var_B = 1, lambda = 1.0001)) |> setDT()
ols_pred_2 <- predict_eqa(data = test_data_cs_2, new_data = test_data_eq_2, method = "ols", level = 0.95, imprecision_estimates = list(Var_A = 3, Var_B = 1, lambda = 1.0001)) |> setDT()
ols_pred_3 <- predict_eqa(data = test_data_cs_3, new_data = test_data_eq_3, method = "ols", level = 0.95, imprecision_estimates = list(Var_A = 3, Var_B = 1, lambda = 1.0001)) |> setDT()

test_that("Check if OLS prediction intervals are wider", {
  expect_true(object = all(ols_pred_1$upr - ols_pred_1$lwr > ss_pred_1$upr - ss_pred_1$lwr))
  expect_true(object = all(ols_pred_2$upr - ols_pred_2$lwr > ss_pred_2$upr - ss_pred_2$lwr))
  expect_true(object = all(ols_pred_3$upr - ols_pred_3$lwr > ss_pred_3$upr - ss_pred_3$lwr))
})

# Check if order is kept

# Draw from sample not necessarily in increasing order
xs <- sample(x = c(1.5, 2:10, 10.5, 11), size = 12, replace = F)

# Update eq_parameters
eq_parameters$obs_tau <- xs

# Simulate evaluated material data
test_data_eq_4 <- sim_eqa_data(parameters = eq_parameters, type = 3) |> setDT()
test_data_eq_4[, SampleID := paste0("EQAM_", SampleID)]

# Prediction
ss_pred_4 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_4, df = 5, level = 0.95, rounding = 6)

test_that("Check ordering of output without NA values", {
  expect_true(object = all(order(test_data_eq_4$MP_B) == order(ss_pred_4$MP_B)))
})

# Assign three of the values as NA values
test_data_eq_5 <- test_data_eq_4
na_ids <- sample(1:12, size = 3, replace = FALSE)
test_data_eq_5$MP_B[na_ids] <- NA

# Prediction
ss_pred_5 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_5)

test_that("Check ordering of output with NA values", {
  expect_true(object = all(order(test_data_eq_5$MP_A) == order(ss_pred_5$MP_A)))
})

# Update eq_parameters
eq_parameters$obs_tau <- seq(from = 2, to = 10, length.out = 1e3)

# Simulate evaluated material data
test_data_eq_6 <- sim_eqa_data(parameters = eq_parameters, type = 2)

# Remove SampleID
test_data_eq_6$SampleID <- NULL

test_that("Check if runs if SampleID missing",{
  expect_no_error(predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_6, df = 5, level = 0.95))
  expect_no_warning(predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_6, df = 5, level = 0.95))
})

ss_pred_6 <- predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_6, df = 5, level = 0.95, rounding = 6)
ss_pred_6_ref <- predict(ss_fit_2_ref, x = test_data_eq_6$MP_B)$y

test_that("Check stuff when SampleID missing",{
  expect_true(object = all(order(test_data_eq_6$MP_B) == order(ss_pred_6$MP_B)))
  expect_true(!is.null(test_data_eq_6$MP_B))
  expect_true(all(as.integer(ss_pred_6$MP_A > ss_pred_6$lwr & ss_pred_6$MP_A < ss_pred_6$upr) == ss_pred_6$inside))
  expect_true(all(abs(ss_pred_6$prediction - ss_pred_6_ref) < 0.01))
})
















