library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

# Testing on "real data"
real_data <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

# Test 1.1.1 - 1.4.4
test_that(desc = "Check structure of smoothing_spline_diagnostics(data, only_optimal_dfs = TRUE)", code = {

  # Setting seed
  set.seed(99)

  # Test data
  test_data_1 <- real_data[comparison == unique(real_data$comparison)[[1]]]
  test_data_2 <- real_data[comparison == unique(real_data$comparison)[[3]]]
  test_data_3 <- real_data[comparison == unique(real_data$comparison)[[9]]]

  # Actual output
  actual_1 <- smoothing_spline_diagnostics(data = test_data_1, only_optimal_dfs = TRUE, weighted = FALSE, exclude_not_restricted = TRUE)
  actual_2 <- smoothing_spline_diagnostics(data = test_data_2, only_optimal_dfs = TRUE, weighted = TRUE, exclude_not_restricted = TRUE)
  actual_3 <- smoothing_spline_diagnostics(data = test_data_3, only_optimal_dfs = TRUE, weighted = FALSE, exclude_not_restricted = FALSE)
  actual_4 <- smoothing_spline_diagnostics(data = test_data_2, only_optimal_dfs = TRUE, weighted = TRUE, exclude_not_restricted = FALSE)

  # Check if data.table class
  expect_true(object = is.data.table(actual_1),
              label = "[TEST 1.1.1]")
  expect_true(object = is.data.table(actual_2),
              label = "[TEST 1.1.2]")
  expect_true(object = is.data.table(actual_3),
              label = "[TEST 1.1.3]")
  expect_true(object = is.data.table(actual_4),
              label = "[TEST 1.1.4]")

  # Check number of rows and columns
  expect_true(object = nrow(actual_1) == 3L & ncol(actual_1) == 4L,
              label = "[TEST 1.2.1]")
  expect_true(object = nrow(actual_2) == 6L & ncol(actual_2) == 5L,
              label = "[TEST 1.2.2]")
  expect_true(object = nrow(actual_3) == 6L & ncol(actual_3) == 5L,
              label = "[TEST 1.2.3]")
  expect_true(object = nrow(actual_4) == 12L & ncol(actual_4) == 6L,
              label = "[TEST 1.2.4]")

  # Check column names
  expect_named(object = actual_1,
               expected = c("df", "method", "smudged", "cv"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.3.1]")
  expect_named(object = actual_2,
               expected = c("df", "method", "weighted", "smudged", "cv"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.3.2]")
  expect_named(object = actual_3,
               expected = c("df", "method", "restricted", "smudged", "cv"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.3.3]")
  expect_named(object = actual_4,
               expected = c("df", "method", "restricted", "weighted", "smudged", "cv"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.3.4]")

  # Check method column contents
  expect_true(object = sum(actual_1$method %in% c("LOOCV", "GCV")) == 3,
              label = "[TEST 1.4.1]")
  expect_true(object = sum(actual_2$method %in% c("LOOCV", "GCV")) == 6,
              label = "[TEST 1.4.2]")
  expect_true(object = sum(actual_3$method %in% c("LOOCV", "GCV")) == 6,
              label = "[TEST 1.4.3]")
  expect_true(object = sum(actual_4$method %in% c("LOOCV", "GCV")) == 12,
              label = "[TEST 1.4.4]")

})

# Test 2.1.1 -
# Include tests on numeric values in smoothing_spline_diagnostics(data, only_optimal_dfs = TRUE)

# Test 3.1.1 -
# Include tests on structure of smoothing_spline_diagnostics(data, only_optimal_dfs = FALSE)






