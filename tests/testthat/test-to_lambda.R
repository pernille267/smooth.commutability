library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)
library(splines)

test_data <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0), type = 2)

x_orig <- test_data$MP_B[order(test_data$MP_B)]
x_unit <- (x_orig - min(x_orig)) / diff(range(x_orig))
k_inte <- calculate_interior_knots(x_unit)
B <- spline.des(knots = c(rep(0, 4L), k_inte, rep(1, 4L)), x = x_unit, ord = 4, derivs = 0, outer.ok = TRUE)$design
BTB <- crossprod(B, B)
Omega <- penalty_matrix(x_min = 0, x_max = 1, interior_knots = k_inte)

true_lambda_1 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 3, all.knots = c(0, k_inte, 1))$lambda
true_lambda_2 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 4, all.knots = c(0, k_inte, 1))$lambda
true_lambda_3 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 5, all.knots = c(0, k_inte, 1))$lambda
true_lambda_4 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 6, all.knots = c(0, k_inte, 1))$lambda
true_lambda_5 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 7.5, all.knots = c(0, k_inte, 1))$lambda
true_lambda_6 <- smooth.spline(x = test_data$MP_B, y = test_data$MP_A, df = 13, all.knots = c(0, k_inte, 1))$lambda

actual_lambda_1 <- to_lambda(df = 3, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)
actual_lambda_2 <- to_lambda(df = 4, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)
actual_lambda_3 <- to_lambda(df = 5, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)
actual_lambda_4 <- to_lambda(df = 6, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)
actual_lambda_5 <- to_lambda(df = 7.5, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)
actual_lambda_6 <- to_lambda(df = 13, weights = rep(1, 25), B = B, BTB = BTB, Omega = Omega)

test_that("Lambda corresponds with that of smooth.spline with a relative error margin smaller than 0.5%", {
  expect_true(abs(actual_lambda_1 - true_lambda_1) / true_lambda_1 * 100 < 0.5)
  expect_true(abs(actual_lambda_2 - true_lambda_2) / true_lambda_2 * 100 < 0.5)
  expect_true(abs(actual_lambda_3 - true_lambda_3) / true_lambda_3 * 100 < 0.5)
  expect_true(abs(actual_lambda_4 - true_lambda_4) / true_lambda_4 * 100 < 0.5)
  expect_true(abs(actual_lambda_5 - true_lambda_5) / true_lambda_5 * 100 < 0.5)
  expect_true(abs(actual_lambda_6 - true_lambda_6) / true_lambda_6 * 100 < 0.5)
})
