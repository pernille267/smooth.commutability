library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

f2 <- function(x){
  return(x + 0.05 * exp(0.16 * x^1.35))
}

inverse <- function(f, lower = 0, upper = 1e3){
  function(x) {
    uniroot(function(y) f(y) - x, lower = lower, upper = upper)$root
  }
}

f2_inv <- inverse(f2)

derivate <- function(f, x){
  h <- .Machine$double.eps^(1/(1+2))
  return((do.call(what = f, args = list(x = x + h)) - do.call(what = f, args = list(x = x - h))) / 2 / h)
}

df2_orig <- function(x, squared = FALSE){
  return(derivate(f = f2, x = x)^(1 + as.integer(squared)))
}

df2_inv <- function(x, squared = FALSE){
  return(derivate(f = f2_inv, x = x)^(1 + as.integer(squared)))
}

df2_orig_mean <- function(lower = 2, upper = 10, squared = FALSE){
  cfun <- function(a){
    df2_orig(x = a, squared = squared) * 1/(upper - lower)
  }
  return(integrate(cfun, lower = lower, upper = upper)$value)
}

df2_inv_mean <- function(lower = 2, upper = 10, squared = FALSE){
  x <- runif(n = 1e6, min = lower, max = upper)
  dens <- density(sapply(X = x, FUN = f2), from = f2(lower) + 0.45, to = f2(upper) - 0.45)
  ss_obj <- smooth.spline(x = dens$x, y = dens$y)
  cfun <- function(a){
    sapply(a, function(b) df2_inv(x = b, squared = squared) * predict(ss_obj, x = b)$y)
  }
  dfun <- function(a){
    sapply(a, function(b) predict(ss_obj, x = b)$y)
  }
  #print(integrate(dfun, lower = f2(lower), upper = f2(upper))$value)
  return(integrate(cfun, lower = f2(lower), upper = f2(upper))$value)
}

true_inv_slope_mean <- df2_inv_mean(squared = FALSE, lower = 2, upper = 10)
true_inv_squared_slope_mean <- df2_inv_mean(squared = TRUE, lower = 2, upper = 10)
true_slope_mean <- df2_orig_mean(squared = FALSE, lower = 2, upper = 10)
true_squared_slope_mean <- df2_orig_mean(squared = TRUE, lower = 2, upper = 10)


pars_1 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0)
pars_2 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, qpos = 1, qran = 0.5, mmax = 10)
pars_3 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, b0 = 0, b1 = 1.05)
pars_4 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.05, cvy = 0.01, cve = 0)
test_data_1 <- simulate_eqa_data2(parameters = pars_1, type = 2, AR = TRUE) |> setDT()
test_data_2 <- simulate_eqa_data(parameters = pars_2) |> setDT()
test_data_3 <- simulate_eqa_data(parameters = pars_3) |> setDT()
test_data_4 <- simulate_eqa_data(parameters = pars_4) |> setDT()

test_that("Does estimate_zeta_ss correspond with estimate_zeta if df <= 2", {
  expect_equal(estimate_zeta_ss(data = test_data_1, df = 2)$zeta, estimate_zeta(test_data_1)$zeta)
  expect_equal(estimate_zeta_ss(data = test_data_2, df = 1)$zeta, estimate_zeta(test_data_2)$zeta)
  expect_equal(estimate_zeta_ss(data = test_data_3, df = 1.5)$zeta, estimate_zeta(test_data_3)$zeta)
})

test_that("Estimated zeta values are positive", {
  expect_gt(object = estimate_zeta_ss(data = test_data_1, df = runif(1, 3, 10))$zeta, expected = 0)
  expect_gt(object = estimate_zeta_ss(data = test_data_2, df = runif(1, 3, 10))$zeta, expected = 0)
  expect_gt(object = estimate_zeta_ss(data = test_data_3, df = runif(1, 3, 10))$zeta, expected = 0)
})

test_that("Increasing df from 2 decrease zeta if there is a non-linear relationship", {
  expect_lt(object = estimate_zeta_ss(data = test_data_1, df = 3)$zeta, expected = estimate_zeta_ss(data = test_data_1, df = 2)$zeta)
  expect_lt(object = estimate_zeta_ss(data = test_data_2, df = 3)$zeta, expected = estimate_zeta_ss(data = test_data_2, df = 2)$zeta)
})

many_zetas <- replicate(n = 1e3, expr = estimate_zeta_ss(data = simulate_eqa_data2(parameters = pars_4, type = 2, AR = TRUE, shift = FALSE), df = 5, simple_output = FALSE), simplify = FALSE) |> rbindlist()
emp_mean_slopes <- many_zetas[many_zetas$shift]$slopes |> mean()
emp_mean_squared_slopes <- many_zetas[many_zetas$shift]$sq_slopes |> mean()

test_that("Test inverted slopes", {
  expect_true(abs(emp_mean_slopes - true_inv_slope_mean) / true_inv_slope_mean * 100 < 2.5, label = "[Average slope match theoretical slope]")
  expect_true(abs(emp_mean_squared_slopes - true_inv_squared_slope_mean) / true_inv_squared_slope_mean * 100 < 5, label = "[Average squared slope match theoretical sqaured slope]")
})

many_zetas <- replicate(n = 1e3, expr = estimate_zeta_ss(data = simulate_eqa_data2(parameters = pars_1, type = 2, AR = TRUE, shift = FALSE), df = 5, simple_output = FALSE), simplify = FALSE) |> rbindlist()
emp_mean_slopes <- many_zetas[!many_zetas$shift]$slopes |> mean()
emp_mean_squared_slopes <- many_zetas[!many_zetas$shift]$sq_slopes |> mean()

test_that("Test original slopes", {
  expect_true(abs(emp_mean_slopes - true_slope_mean) / true_slope_mean * 100 < 2.5, label = "[Average slope match theoretical slope]")
  expect_true(abs(emp_mean_squared_slopes - true_squared_slope_mean) / true_squared_slope_mean * 100 < 5, label = "[Average squared slope match theoretical sqaured slope]")
})

test_that("Test simple output", {
  expect_length(estimate_zeta_ss(data = test_data_1, df = 5, simple_output = TRUE), n = 1)
  expect_length(estimate_zeta_ss(data = test_data_1, df = 5, simple_output = FALSE), n = 8)
  expect_named(estimate_zeta_ss(data = test_data_1, df = 5, simple_output = TRUE), expected = "zeta")
  expect_named(object = estimate_zeta_ss(data = test_data_1, df = 5, simple_output = FALSE),
               expected = c("zeta", "irr_var", "tot_var", "var_h", "var_v", "slopes", "sq_slopes", "shift"))
})

test_that("Test that irr_var is correctly calculated", {
  expected_irr_var <- many_zetas[!many_zetas$shift]$var_v + many_zetas[!many_zetas$shift]$sq_slopes * many_zetas[!many_zetas$shift]$var_h
  actual_irr_var <- many_zetas[!many_zetas$shift]$irr_var
  true_irr_var <- (0.01 * 6)^2 + true_squared_slope_mean * (0.01 * 6)^2
  expect_true(all(expected_irr_var == actual_irr_var))
  expect_true(abs(mean(actual_irr_var) - true_irr_var) / true_irr_var * 100 < 5)
})





