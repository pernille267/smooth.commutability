# Repr.
set.seed(99)

# Test data
test_data_1 <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10), type = 2, AR = TRUE) |> setDT()
test_data_2 <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

# Test if global_precision_estimates2() works as intended
test_that("Check global_precision_estimates2()", code = {

  # Actual and expected output
  actual_1 <- as.data.table(global_precision_estimates2(test_data_1))
  actual_2 <- test_data_2[, global_precision_estimates2(.SD), by = comparison]
  expected_1 <- test_data_1[, list(Var_A = var(MP_A, na.rm = TRUE),
                                   Var_B = var(MP_B, na.rm = TRUE),
                                   m_A = mean(MP_A, na.rm = TRUE),
                                   m_B = mean(MP_B, na.rm = TRUE)),
                            by = SampleID]
  expected_1 <- expected_1[, list(Var_A = mean(Var_A, na.rm = TRUE),
                                  Var_B = mean(Var_B, na.rm = TRUE),
                                  CV_A = sqrt(mean(Var_A, na.rm = TRUE)) / mean(m_A, na.rm = TRUE),
                                  CV_B = sqrt(mean(Var_B, na.rm = TRUE)) / mean(m_B, na.rm = TRUE),
                                  lambda = mean(Var_A, na.rm = TRUE) / mean(Var_B, na.rm = TRUE))]
  expected_2 <- test_data_2[, list(Var_A = var(MP_A, na.rm = TRUE),
                                   Var_B = var(MP_B, na.rm = TRUE),
                                   m_A = mean(MP_A, na.rm = TRUE),
                                   m_B = mean(MP_B, na.rm = TRUE)),
                            by = list(comparison, SampleID)]
  expected_2 <- expected_2[, list(Var_A = mean(Var_A, na.rm = TRUE),
                                  Var_B = mean(Var_B, na.rm = TRUE),
                                  CV_A = sqrt(mean(Var_A, na.rm = TRUE)) / mean(m_A, na.rm = TRUE),
                                  CV_B = sqrt(mean(Var_B, na.rm = TRUE)) / mean(m_B, na.rm = TRUE),
                                  lambda = mean(Var_A, na.rm = TRUE) / mean(Var_B, na.rm = TRUE)),
                           by = comparison]

  # Check if equal
  expect_equal(object = actual_1, expected = expected_1)
  expect_equal(object = actual_2, expected = expected_2)

})

# Test if resample_samples2() works as intended
test_that("Check resample_samples2()", code = {

  # Relevant components
  test_data_1_mor <- test_data_1[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = SampleID]
  test_data_1_vor <- test_data_1[, list(MP_A = var(MP_A), MP_B = var(MP_B)), by = SampleID]

  # Actual output
  actual_1 <- as.data.table(resample_samples2(test_data_1))
  actual_1$SampleID <- as.integer(actual_1$SampleID)
  actual_1$ReplicateID <- as.integer(actual_1$Replicate)
  actual_1_mor <- actual_1[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = SampleID]
  actual_1_vor <- actual_1[, list(MP_A = var(MP_A), MP_B = var(MP_B)), by = SampleID]

  # Test if resampled data exists in original data
  expect_true(object = all(actual_1$MP_A %in% test_data_1$MP_A))
  expect_true(object = all(actual_1$MP_B %in% test_data_1$MP_B))

  # Test if data is resampled at SampleID-level (*Pseudo-tests)
  expect_true(object = all(actual_1_mor$MP_A %in% test_data_1_mor$MP_A))
  expect_true(object = all(actual_1_mor$MP_B %in% test_data_1_mor$MP_B))
  expect_true(object = all(actual_1_vor$MP_A %in% test_data_1_vor$MP_A))
  expect_true(object = all(actual_1_vor$MP_B %in% test_data_1_vor$MP_B))

  # Test if structure is correct
  expect_named(object = actual_1, expected = names(test_data_1))
  expect_equal(object = unique(actual_1$SampleID), expected = unique(test_data_1$SampleID))
  expect_true(object = all(actual_1$SampleID == test_data_1$SampleID))
})

# Test if resample_fun_of_samples() works as intended
test_that("Check resample_fun_of_samples()", code = {

  # Test data used in this function
  test_data_1_mor <- test_data_1[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = SampleID]
  test_data_1_vor <- test_data_1[, list(MP_A = var(MP_A), MP_B = var(MP_B)), by = SampleID]

  # Actual output
  actual_1 <- as.data.table(resample_fun_of_samples(test_data_1_mor))
  actual_1$SampleID <- as.integer(actual_1$SampleID)
  actual_2 <- as.data.table(resample_fun_of_samples(test_data_1_vor))
  actual_2$SampleID <- as.integer(actual_2$SampleID)

  # Test if resampled data exists in original data
  expect_true(object = all(actual_1$MP_A %in% test_data_1_mor$MP_A))
  expect_true(object = all(actual_1$MP_B %in% test_data_1_mor$MP_B))
  expect_true(object = all(actual_2$MP_A %in% test_data_1_vor$MP_A))
  expect_true(object = all(actual_2$MP_B %in% test_data_1_vor$MP_B))

  # Test if structure is correct
  expect_named(object = actual_1, expected = names(test_data_1_mor))
  expect_equal(object = unique(actual_1$SampleID), expected = unique(test_data_1_mor$SampleID))
  expect_true(object = all(actual_1$SampleID == test_data_1_mor$SampleID))
  expect_named(object = actual_2, expected = names(test_data_1_vor))
  expect_equal(object = unique(actual_2$SampleID), expected = unique(test_data_1_vor$SampleID))
  expect_true(object = all(actual_2$SampleID == test_data_1_vor$SampleID))

})

# Test if resample_fun_of_samples_all() works as intended
test_that("Check resample_fun_of_samples_all()", code = {

  # Test data used in this function
  test_data_2_mor <- test_data_2[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = list(comparison, SampleID)]
  test_data_2_vor <- test_data_2[, list(MP_A = var(MP_A), MP_B = var(MP_B)), by = list(comparison, SampleID)]

  # Actual output
  actual_1 <- as.data.table(resample_fun_of_samples_all(test_data_2_mor))
  actual_1$SampleID <- as.integer(actual_1$SampleID)
  actual_2 <- as.data.table(resample_fun_of_samples_all(test_data_2_vor))
  actual_2$SampleID <- as.integer(actual_2$SampleID)

  # Test if resampled data exists in original data
  for(i in unique(test_data_2_mor$comparison)){
    expect_true(object = all(actual_1[comparison == i]$MP_A %in% test_data_2_mor[comparison == i]$MP_A), label = i)
    expect_true(object = all(actual_1[comparison == i]$MP_B %in% test_data_2_mor[comparison == i]$MP_B), label = i)
    expect_true(object = all(actual_2[comparison == i]$MP_A %in% test_data_2_vor[comparison == i]$MP_A), label = i)
    expect_true(object = all(actual_2[comparison == i]$MP_B %in% test_data_2_vor[comparison == i]$MP_B), label = i)
  }
})

# Test if resample_imprecision() works as intended (*Pseudo-tests)
test_that("Check resample_imprecision()", code = {

  true_sd_CV_A <- 0.1225967
  true_sd_lambda <- 0.3037452

  resampled_CV_A <- replicate(n = 1000, resample_imprecision(test_data_1)$CV_A * 100)
  resampled_lambda <- replicate(n = 1000, resample_imprecision(test_data_1)$lambda)

  expect_true(object = abs(sd(resampled_CV_A, na.rm = TRUE) - true_sd_CV_A) < 0.1)
  expect_true(object = abs(sd(resampled_lambda, na.rm = TRUE) - true_sd_lambda) < 0.2)

})
