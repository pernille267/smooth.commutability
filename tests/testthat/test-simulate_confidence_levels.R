library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

set.seed(99)

parameters_1 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, qpos = 1, qran = 0.3, mmax = 5, df = 5)
parameters_2 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, qpos = 1, qran = 0.3, mmax = 5, df_max = 5)
parameters_3 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0, df = 7.5, type = 2)

sim_conf_lvl_1 <- simulate_confidence_level(N = 50, n = 2, m = 1, parameters = parameters_1, level = 0.95, calculate_sem = TRUE)
sim_conf_lvl_2 <- simulate_confidence_level(N = 100, n = 1, m = 1, parameters = parameters_2, level = 0.95, calculate_sem = TRUE)
sim_conf_lvl_3 <- simulate_confidence_level(N = 500, n = 1, m = 1, parameters = parameters_3, level = 0.95, calculate_sem = TRUE)

true_sem_conf_lvl_3 <- 0.01008897
true_conf_lvl_3 <- 0.9422844

test_that("Checks if confidence levels are in valid range", {
  expect_true(object = sim_conf_lvl_2$confidence_level >= 0 & sim_conf_lvl_2$confidence_level <= 1)
  expect_true(object = sim_conf_lvl_3$confidence_level >= 0 & sim_conf_lvl_3$confidence_level <= 1)
})

test_that("Checks if sem is correctly calculated", {
  expect_true(abs(sqrt(sim_conf_lvl_2$confidence_level * (1 - sim_conf_lvl_2$confidence_level) / sim_conf_lvl_2$number_of_observations) - sim_conf_lvl_2$sem) < 1e-3)
  expect_true(abs(sqrt(sim_conf_lvl_3$confidence_level * (1 - sim_conf_lvl_3$confidence_level) / sim_conf_lvl_3$number_of_observations) - sim_conf_lvl_3$sem) < 1e-3)
})

test_that("Checks if sem is approximately correct", {
  expect_true(abs(sim_conf_lvl_3$sem - true_sem_conf_lvl_3) / true_sem_conf_lvl_3 < 0.25)
})

test_that("Checks if number_of_observations are correct", {
  expect_true(100 - sim_conf_lvl_2$number_of_observations >= 0 & 100 - sim_conf_lvl_2$number_of_observations <= 2)
  expect_true(500 - sim_conf_lvl_3$number_of_observations >= 0 & 500 - sim_conf_lvl_3$number_of_observations <= 10)
})

test_that("Checks replication of simulatd confidence levels", {
  expect_true(nrow(sim_conf_lvl_1) == 2)
  expect_true(all(sim_conf_lvl_1$replication %in% 1:2))
})

test_that("Error handling", {
  expect_error(simulate_confidence_level(N = 1, n = 1, m = 1, parameters = parameters_3, level = 0.95), regexp = "N < 2, but")
  expect_error(simulate_confidence_level(N = 1e7L, n = 1, m = 1, parameters = parameters_3, level = 0.95), regexp = "but this will take forever to run")
  expect_error(simulate_confidence_level(N = NA, n = 1, m = 1, parameters = parameters_3, level = 0.95), regexp = "N is not numeric or is missing")
  expect_error(simulate_confidence_level(N = "text?", n = 1, m = 1, parameters = parameters_3, level = 0.95), regexp = "N is not numeric or is missing")
  expect_error(simulate_confidence_level(N = 1e3L, n = 1.5, m = 1, parameters = parameters_3, level = 0.95), regexp = "n is not an integer")
  expect_error(simulate_confidence_level(N = 1e3L, n = 1e6L, m = 1, parameters = parameters_3, level = 0.95), regexp = "but this will take forever to run")
  expect_error(simulate_confidence_level(N = 1e3L, n = 1, m = 1e4L, parameters = parameters_3, level = 0.95), regexp = "but this choice is insane")
  expect_error(simulate_confidence_level(N = 1e3L, n = 1, m = 1, parameters = parameters_3, level = 1.05), regexp = "level is outside the real-open interval")
  expect_error(simulate_confidence_level(N = 1e3L, n = 1, m = 1, parameters = parameters_3, level = 0.95, calculate_sem = "yes please..."), regexp = "calculate_sem is not a non-missing logical value.")
})

parameters_4 <- CJ(n = c(25, 40), R = c(3, 4), cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0, type = 2, df = 7.5)

output_percent <- simulate_confidence_levels(N = 1e2, n = 1, m = 1, parameters = parameters_4, level = 0.95, parallel = FALSE, calculate_sem = TRUE, percent = TRUE)
output_not_percent <- simulate_confidence_levels(N = 1e2, n = 1, m = 1, parameters = parameters_4, level = 0.95, parallel = FALSE, calculate_sem = TRUE, percent = FALSE)

test_that("Test structure of outputs", {
  expect_true(nrow(output_percent) == 4 & ncol(output_percent) == 3 + ncol(parameters_4))
  expect_true(nrow(output_not_percent) == 4 & ncol(output_not_percent) == 3 + ncol(parameters_4))
  expect_named(output_percent, expected = c(names(parameters_4), "confidence_level", "number_of_observations", "sem"))
  expect_named(output_not_percent, expected = c(names(parameters_4), "confidence_level", "number_of_observations", "sem"))
})

test_that("Test percentages", {
  expect_true(all(output_percent$cvx == output_not_percent$cvx * 100))
  expect_true(all(output_percent$cvy == output_not_percent$cvy * 100))
  expect_true(all(abs(output_percent$confidence_level - output_not_percent$confidence_level * 100) <= 15))
  expect_true(all(abs(output_percent$sem - output_not_percent$sem * 100) <= 3))
})


