library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

simulate_eq_data <- function(x = 5, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10, type = 2, beta0 = 0, beta1 = 1){
  mu_tau <- mean(c(cil, ciu))
  varh <- (mu_tau * cvx)^2 / R
  varv <- (mu_tau * cvy)^2 / R
  m <- length(x)
  h <- rnorm(n = m, mean = 0, sd = sqrt(varh))
  v <- rnorm(n = m, mean = 0, sd = sqrt(varv))
  if(type == 1){
    ny <- x + 0.90 * sin(0.40 * x^1.06) + v
    nx <- x + h
  }
  else if(type == 2){
    ny <- x + 0.05 * exp(0.16 * x^1.35) + v
    nx <- x + h
  }
  else if(type == 3){
    ny <- x - exp(-0.5 * (x - 1.5)^2 / 4) + v
    nx <- x + h
  }
  else{
    ny <- beta0 + beta1 * x + v
    nx <- x + h
  }
  sampleid <- paste0("EQA_", 1:m)
  return(data.table(SampleID = sampleid,
                    MP_A = ny,
                    MP_B = nx))
}


set.seed(99)

test_data_cs_1 <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cve = 0, cil = 2, ciu = 10), type = 1, AR = FALSE) |> setDT()
test_data_cs_2 <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cve = 0, cil = 2, ciu = 10), type = 2, AR = FALSE) |> setDT()
test_data_cs_3 <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cve = 0, cil = 2, ciu = 10), type = 3, AR = FALSE) |> setDT()
test_data_eq_1 <- simulate_eq_data(x = c(1.5, 7, 11), type = 1)
test_data_eq_2 <- simulate_eq_data(x = c(1.5, 7, 11), type = 2)
test_data_eq_3 <- simulate_eq_data(x = c(1.5, 7, 11), type = 3)

ss_fit_1 <- smoothing_spline(data = test_data_cs_1, df = 5)
ss_fit_2 <- smoothing_spline(data = test_data_cs_2, df = 5)
ss_fit_3 <- smoothing_spline(data = test_data_cs_3, df = 5)

ss_fit_1_ref <- smooth.spline(x = test_data_cs_1$MP_B, test_data_cs_1$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_1$interior_knots, 1))
ss_fit_2_ref <- smooth.spline(x = test_data_cs_2$MP_B, test_data_cs_2$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_2$interior_knots, 1))
ss_fit_3_ref <- smooth.spline(x = test_data_cs_3$MP_B, test_data_cs_3$MP_A, df = 5, keep.stuff = TRUE, all.knots = c(0, ss_fit_3$interior_knots, 1))

ss_pred_1 <- predict_smoothing_spline(data = test_data_cs_1, new_data = test_data_eq_1, df = 5, level = 0.95, rounding = 6)
ss_pred_2 <- predict_smoothing_spline(data = test_data_cs_2, new_data = test_data_eq_2, df = 5, level = 0.95, rounding = 6)
ss_pred_3 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_3, df = 5, level = 0.95, rounding = 6)

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

xs <- sample(x = c(1.5, 2:10, 10.5, 11), size = 12, replace = F)
test_data_eq_4 <- simulate_eq_data(x = xs, type = 3)
ss_pred_4 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_4, df = 5, level = 0.95, rounding = 6)

test_that("Check ordering of output without NA values", {
  expect_true(object = all(order(test_data_eq_4$MP_B) == order(ss_pred_4$MP_B)))
})

test_data_eq_5 <- test_data_eq_4
na_ids <- sample(1:12, size = 3, replace = FALSE)
test_data_eq_5$MP_B[na_ids] <- NA

ss_pred_5 <- predict_smoothing_spline(data = test_data_cs_3, new_data = test_data_eq_5)

test_that("Check ordering of output with NA values", {
  expect_true(object = all(order(test_data_eq_5$MP_A) == order(ss_pred_5$MP_A)))
})

test_data_eq_6 <- simulate_eq_data(x = seq(from = 2, to = 10, length.out = 1e3), type = 2)
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

