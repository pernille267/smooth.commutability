library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

set.seed(99)

parameters_1 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, df = 3)
parameters_2 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, df_max = 3)
parameters_3 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01)

test_that("Check if simulate_inside gives appropriate results", {
  expect_true(simulate_inside(parameters = parameters_1, m = 1, level = 0.95) %in% c(0L, 1L, NA_integer_))
  expect_true(simulate_inside(parameters = parameters_2, m = 1, level = 0.95) %in% c(0L, 1L, NA_integer_))
  expect_error(simulate_inside(parameters = parameters_3, m = 1, level = 0.95))
})

test_that("Check if confidence levels are approximately ok", {
  expect_true(abs(mean(replicate(n = 500, simulate_inside(parameters = parameters_1, m = 1, level = 0.95)), na.rm = TRUE) - 0.95) < 0.05)
  expect_true(abs(mean(replicate(n = 500, simulate_inside(parameters = parameters_1, m = 2, level = 0.95)), na.rm = TRUE) - 0.95) < 0.05)
})

test_sis_1 <- simulate_insides(N = 10, m = 1L, parameters = parameters_1, level = 0.95, attach = FALSE, replication_id = FALSE)
test_sis_2 <- simulate_insides(N = 10, m = 1L, parameters = parameters_1, level = 0.95, attach = TRUE, replication_id = FALSE)
test_sis_3 <- simulate_insides(N = 10, m = 2L, parameters = parameters_1, level = 0.95, attach = TRUE, replication_id = TRUE)
test_sis_4 <- simulate_insides(N = 10, m = 4L, parameters = parameters_1, level = 0.95, attach = FALSE, replication_id = TRUE)

test_that("Check names of 'simulate_inside' outputs", {
  expect_named(test_sis_1, expected = "inside", ignore.order = FALSE, ignore.case = FALSE)
  expect_named(test_sis_2, expected = c("n", "R", "cil", "ciu", "cvx", "cvy", "df", "inside"), ignore.order = FALSE, ignore.case = FALSE)
  expect_named(test_sis_3, expected = c("n", "R", "cil", "ciu", "cvx", "cvy", "df", "replication", "id", "SampleID", "inside"), ignore.order = FALSE, ignore.case = FALSE)
  expect_named(test_sis_4, expected = c("replication", "id", "SampleID", "inside"), ignore.order = FALSE, ignore.case = FALSE)
})

test_that("Check dimensions of 'simulate_inside' outputs", {
  expect_true(nrow(test_sis_1) == 10 & ncol(test_sis_1) == 1)
  expect_true(nrow(test_sis_2) == 10 & ncol(test_sis_2) == 8)
  expect_true(nrow(test_sis_3) == 20 & ncol(test_sis_3) == 11)
  expect_true(nrow(test_sis_4) == 40 & ncol(test_sis_4) == 4)
})

test_that("Check cell values of 'simulate_inside' outputs", {
  expect_true(is.integer(test_sis_1$inside))
  expect_true(is.integer(test_sis_2$inside))
  expect_true(is.integer(test_sis_3$inside))
  expect_true(all(test_sis_3$SampleID %in% 1:2))
  expect_true(all(test_sis_3$replication %in% c(1)))
  expect_true(all(test_sis_3$id %in% 1:10))
  expect_true(is.integer(test_sis_4$inside))
  expect_true(all(test_sis_4$SampleID %in% 1:4))
  expect_true(all(test_sis_4$replication %in% c(1)))
  expect_true(all(test_sis_4$id %in% 1:10))
})

test_that("Check if confidence levels are approximately ok again", {
  expect_true(abs(mean(simulate_insides(N = 500L, m = 1, parameters = parameters_1, level = 0.95)$inside, na.rm = TRUE) - 0.95) < 0.05)
  expect_true(abs(mean(simulate_insides(N = 500L, m = 1, parameters = parameters_2, level = 0.95)$inside, na.rm = TRUE) - 0.95) < 0.05)
})


