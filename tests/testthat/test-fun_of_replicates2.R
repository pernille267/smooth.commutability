# Repr.
set.seed(99)

# Test data
test_data_1 <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10), type = 2, AR = TRUE) |> setDT()
test_data_2 <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

# Test if fun_of_replicates2() works as intended
test_that(desc = "Check fun_of_replicates2()", code = {

  # Check mean of replicates:

  # Actual output
  actual_1 <- as.data.table(fun_of_replicates2(test_data_1))
  actual_1[, SampleID := as.integer(SampleID)]
  actual_2 <- test_data_2[, fun_of_replicates2(.SD), by = comparison]
  actual_2[, SampleID := as.integer(SampleID)]

  # Expected output
  expected_1 <- test_data_1[, list(MP_A = mean(MP_A, na.rm = TRUE),
                                   MP_B = mean(MP_B, na.rm = TRUE)),
                            by = SampleID]
  expected_2 <- test_data_2[, list(MP_A = mean(MP_A),
                                   MP_B = mean(MP_B)),
                            by = list(comparison, SampleID)]

  # Check if equal
  expect_equal(object = actual_1, expected = expected_1)
  expect_equal(object = actual_2, expected = expected_2)

  # Check variance of replicates:

  # Actual output
  actual_1 <- as.data.table(fun_of_replicates2(test_data_1, "var"))
  actual_1[, SampleID := as.integer(SampleID)]
  actual_2 <- test_data_2[, fun_of_replicates2(.SD, "var"), by = comparison]
  actual_2[, SampleID := as.integer(SampleID)]

  # Expected output
  expected_1 <- test_data_1[, list(MP_A = var(MP_A, na.rm = TRUE),
                                   MP_B = var(MP_B, na.rm = TRUE)),
                            by = SampleID]
  expected_2 <- test_data_2[, list(MP_A = var(MP_A, na.rm = TRUE),
                                   MP_B = var(MP_B, na.rm = TRUE)),
                            by = list(comparison, SampleID)]

  # Check if equal
  expect_equal(object = actual_1, expected = expected_1)
  expect_equal(object = actual_2, expected = expected_2)

  # Check standard deviation of replicates:

  # Actual output
  actual_1 <- as.data.table(fun_of_replicates2(test_data_1, "sd"))
  actual_1[, SampleID := as.integer(SampleID)]
  actual_2 <- test_data_2[, fun_of_replicates2(.SD, "sd"), by = comparison]
  actual_2[, SampleID := as.integer(SampleID)]

  # Expected output
  expected_1 <- test_data_1[, list(MP_A = sd(MP_A, na.rm = TRUE),
                                   MP_B = sd(MP_B, na.rm = TRUE)),
                            by = SampleID]
  expected_2 <- test_data_2[, list(MP_A = sd(MP_A, na.rm = TRUE),
                                   MP_B = sd(MP_B, na.rm = TRUE)),
                            by = list(comparison, SampleID)]

  # Check if equal
  expect_equal(object = actual_1, expected = expected_1)
  expect_equal(object = actual_2, expected = expected_2)

  # Check coefficient of variation of replicates:

  # Actual output
  actual_1 <- as.data.table(fun_of_replicates2(test_data_1, "cv"))
  actual_1[, SampleID := as.integer(SampleID)]
  actual_2 <- test_data_2[, fun_of_replicates2(.SD, "cv"), by = comparison]
  actual_2[, SampleID := as.integer(SampleID)]

  # Expected output
  expected_1 <- test_data_1[, list(MP_A = sd(MP_A, na.rm = TRUE) / mean(MP_A, na.rm = TRUE),
                                   MP_B = sd(MP_B, na.rm = TRUE) / mean(MP_B, na.rm = TRUE)),
                            by = SampleID]
  expected_2 <- test_data_2[, list(MP_A = sd(MP_A, na.rm = TRUE) / mean(MP_A, na.rm = TRUE),
                                   MP_B = sd(MP_B, na.rm = TRUE) / mean(MP_B, na.rm = TRUE)),
                            by = list(comparison, SampleID)]

  # Check if equal
  expect_equal(object = actual_1, expected = expected_1)
  expect_equal(object = actual_2, expected = expected_2)

})


