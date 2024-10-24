library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

test_data_1 <- fread(file = "~/Packages/test data smooth.commutability/fictive_glucose_cs_data.csv")
test_data_2 <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

test_that("fix_predictor_and_response returns same names as input", {
  expect_equal(names(test_data_1), names(fix_predictor_and_response(test_data_1)))
  expect_equal(names(test_data_2), names(fix_predictor_and_response(test_data_2)))
})

test_that("fix_predictor_and_response returns same number of rows as input", {
  expect_equal(nrow(test_data_1), nrow(fix_predictor_and_response(test_data_1)))
  expect_equal(nrow(test_data_2), nrow(fix_predictor_and_response(test_data_2)))
})

test_that("fix_predictor_and_response returns a data.table object", {
  expect_s3_class(fix_predictor_and_response(test_data_1), class = "data.table")
  expect_s3_class(fix_predictor_and_response(test_data_2), class = "data.table")
})

test_that("fix_predictor_and_response returns output with shift column if prompted", {
  expect_true("shift" %in% names(fix_predictor_and_response(test_data_1, NULL, TRUE)))
  expect_true("shift" %in% names(fix_predictor_and_response(test_data_2, NULL, TRUE)))
})

expected_shifts_1 <- test_data_1[, global_precision_estimates(.SD), by = comparison]$lambda < 0.5
expected_shifts_2 <- test_data_2[, global_precision_estimates(.SD), by = comparison]$lambda < 0.5

test_that("fix_predictor_and_response denote shifts correctly", {
  actual_shifts_1 <- fix_predictor_and_response(test_data_1, NULL, TRUE)[, list(shift = shift[1]), by = comparison]$shift
  actual_shifts_2 <- fix_predictor_and_response(test_data_2, NULL, TRUE)[, list(shift = shift[1]), by = comparison]$shift
  expect_equal(actual_shifts_1, expected_shifts_1)
  expect_equal(actual_shifts_2, expected_shifts_2)
})

test_that("fix_predictor_and_response changes names of comparisons if shifted", {
  expected_comparisons_1 <- ifelse(expected_shifts_1, sub("(.*) - (.*)", "\\2 - \\1", unique(test_data_1$comparison)), unique(test_data_1$comparison))
  expected_comparisons_2 <- ifelse(expected_shifts_2, sub("(.*) - (.*)", "\\2 - \\1", unique(test_data_2$comparison)), unique(test_data_2$comparison))
  actual_comparisons_1 <- unique(fix_predictor_and_response(test_data_1)$comparison)
  actual_comparisons_2 <- unique(fix_predictor_and_response(test_data_2)$comparison)
  expect_equal(actual_comparisons_1, expected_comparisons_1)
  expect_equal(actual_comparisons_2, expected_comparisons_2)
})

test_that("fix_predictor_and_response change contents of MP_A and MP_B if shifted", {
  expected_test_data_1_shifted <- unique(test_data_1$comparison)[expected_shifts_1]
  expected_test_data_2_shifted <- unique(test_data_2$comparison)[expected_shifts_2]
  actual_test_data_1_shifted <- unique(fix_predictor_and_response(test_data_1)$comparison)[expected_shifts_1]
  actual_test_data_2_shifted <- unique(fix_predictor_and_response(test_data_2)$comparison)[expected_shifts_2]

  for(i in 1:length(expected_test_data_1_shifted)){
    expected_1 <- test_data_1[comparison == expected_test_data_1_shifted[i],]
    actual_1 <- fix_predictor_and_response(test_data_1)[comparison == actual_test_data_1_shifted[i]]
    expect_equal(expected_1$MP_A, actual_1$MP_B)
    expect_equal(expected_1$MP_B, actual_1$MP_A)
  }

  for(i in 1:length(expected_test_data_2_shifted)){
    expected_2 <- test_data_2[comparison == expected_test_data_2_shifted[i],]
    actual_2 <- fix_predictor_and_response(test_data_2)[comparison == actual_test_data_2_shifted[i]]
    expect_equal(expected_1$MP_A, actual_1$MP_B)
    expect_equal(expected_1$MP_B, actual_1$MP_A)
  }

})

