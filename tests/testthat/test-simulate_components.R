library(fasteqa)
library(data.table)
library(testthat)

# Reproducibility
set.seed(99)

parameters_1 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 1e-20, cvy = 0.01, df_max = 5)
parameters_2 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 1e-20, cvy = 0.01, type = 2, df_max = 5, obs_tau = "RANDOM")
parameters_3 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 1e-20, cvy = 0.01, type = 1, df_max = 5, obs_tau = c(3, 6, 9))
parameters_4 <- CJ(n = c(25, 40), R = c(3, 4), cil = 2, ciu = 10, cvx = 1e-20, cvy = 0.01, type = c(1, 2), df_max = 5)

test_that(desc = "Check simulate_components()", code = {
  # Function output (w. attempt_fast = TRUE)
  actual_1 <- simulate_components(parameters = parameters_1, shift = FALSE, attempt_fast = TRUE)
  actual_2 <- simulate_components(parameters = parameters_2, shift = FALSE, attempt_fast = TRUE)
  actual_3 <- simulate_components(parameters = parameters_3, shift = FALSE, attempt_fast = TRUE)

  # Tests 1.1.1 - 1.1.3
  expect_named(actual_1, expected = c("var_eps", "df"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 1.1.1]")
  expect_named(actual_2, expected = c("var_eps", "var_pred", "var_pred_error", "df", "t"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 1.1.2]")
  expect_equal(length(actual_3), expected = 9, label = "[TEST 1.1.3]")

  # Tests 1.2.1 - 1.2.2
  expect_lte(actual_1[2], parameters_1$df_max + 1e-1, label = "[TEST 1.2.1]")
  expect_lte(actual_2[4], parameters_2$df_max + 1e-1, label = "[TEST 1.2.2]")

  # Function output (w. attempt_fast = FALSE)
  actual_1 <- simulate_components(parameters = parameters_1, shift = FALSE, attempt_fast = FALSE)
  actual_2 <- simulate_components(parameters = parameters_2, shift = FALSE, attempt_fast = FALSE)
  actual_3 <- simulate_components(parameters = parameters_3, shift = FALSE, attempt_fast = FALSE)

  # Tests 1.3.1 - 1.3.3
  expect_named(actual_1, expected = c("var_eps", "df"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 1.3.1]")
  expect_named(actual_2, expected = c("var_eps", "var_pred", "var_pred_error", "df", "t"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 1.3.2]")
  expect_equal(length(actual_3), expected = 9, label = "[TEST 1.3.3]")

  # Tests 1.4.1 - 1.4.2
  expect_lte(actual_1[2], parameters_1$df_max + 1e-1, label = "[TEST 1.4.1]")
  expect_lte(actual_2[4], parameters_2$df_max + 1e-1, label = "[TEST 1.4.2]")
})

test_that(desc = "Check replicate_simulate_components()", code = {
  # Function output (one set of parameters)
  actual_1 <- replicate_simulate_components(parameters = parameters_1, N = 5, parallel = FALSE)
  actual_2 <- replicate_simulate_components(parameters = parameters_2, N = 5, parallel = FALSE)
  actual_3 <- replicate_simulate_components(parameters = parameters_3, N = 5, parallel = FALSE)

  # Tests 2.1.1 - 2.1.3
  expect_named(actual_1, expected = c(names(parameters_1), "m", "x", "component", "value"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 2.1.1]")
  expect_named(actual_2, expected = c(names(parameters_2), "m", "x", "component", "value"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 2.1.2]")
  expect_named(actual_3, expected = c(names(parameters_3), "m", "x", "component", "value"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 2.1.3]")

  # Tests 2.2.1 - 2.2.9
  expect_equal(actual_1$m, rep(1, 5 * 2), label = "[TEST 2.2.1]")
  expect_equal(actual_1$x, rep(1:5, each = 2), label = "[TEST 2.2.2]")
  expect_equal(actual_1$component, rep(c("var_eps", "df"), times = 5), label = "[TEST 2.2.3]")
  expect_equal(actual_2$m, rep(1, 5 * 5), label = "[TEST 2.2.4]")
  expect_equal(actual_2$x, rep(1:5, each = 5), label = "[TEST 2.2.5]")
  expect_equal(actual_2$component, rep(c("var_eps", "var_pred", "var_pred_error", "df", "t"), times = 5), label = "[TEST 2.2.6]")
  expect_equal(actual_3$m, rep(1:3, times = 5 * 3), label = "[TEST 2.2.7]")
  expect_equal(actual_3$x, rep(1:5, each = 3 * 3), label = "[TEST 2.2.8]")
  expect_equal(actual_3$component, rep(rep(c("var_pred", "var_pred_error", "t"), each = 3), times = 5), label = "[TEST 2.2.9]")

  # Function output (several sets of parameters)
  actual_4 <- replicate_simulate_components(parameters = parameters_4, N = 5, parallel = FALSE)

  # Tests 2.3.1 - 2.3.3
  expect_named(actual_4, expected = c(names(parameters_4), "m", "x", "component", "value"), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 2.3.1]")
  expect_true(is.data.table(actual_4), label = "[TEST 2.3.2]")
  expect_true(nrow(actual_4) == nrow(parameters_4) * 2 * 5, label = "[TEST 2.3.3]")
})

test_that(desc = "Check simulate_components_statistics()", code = {
  # Function output (one set of parameters)
  actual_1 <- simulate_components_statistics(parameters = parameters_1, N = 50, parallel = FALSE)
  actual_2 <- simulate_components_statistics(parameters = parameters_2, N = 50, parallel = FALSE)
  actual_3 <- simulate_components_statistics(parameters = parameters_3, N = 50, parallel = FALSE)

  # Tests 3.1.1 - 3.1.3
  expected_stats <- c("min", "q1", "q2.5", "q5", "q10", "Q1", "median", "Q3", "q90", "q95", "q97.5", "q99", "max",
                      "mean", "var", "sd", "mad", "skewness", "kurtosis")
  expect_named(actual_1, expected = c(names(parameters_1), "component", expected_stats), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 3.1.1]")
  expect_named(actual_2, expected = c(names(parameters_2), "component", expected_stats), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 3.1.2]")
  expect_named(actual_3, expected = c(names(parameters_3), "component", expected_stats), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 3.1.3]")

  # Tests 3.2.1 - 3.2.3
  expect_true(nrow(actual_1) == 2, label = "[TEST 3.2.1]")
  expect_true(nrow(actual_2) == 5, label = "[TEST 3.2.2]")
  expect_true(nrow(actual_3) == 9, label = "[TEST 3.2.3]")

  # Function output (several sets of parameters)
  actual_4 <- simulate_components_statistics(parameters = parameters_4, N = 20, parallel = FALSE)

  # Tests 3.3.1 - 3.3.3
  expect_named(actual_4, expected = c(names(parameters_4), "component", expected_stats), ignore.order = FALSE, ignore.case = FALSE, label = "[TEST 3.3.1]")
  expect_true(is.data.table(actual_4), label = "[TEST 3.3.2]")
  expect_true(nrow(actual_4) == nrow(parameters_4) * 2, label = "[TEST 3.3.3]")
})
