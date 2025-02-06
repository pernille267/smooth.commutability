# Repr.
set.seed(99)

# parameters
parameters_1 <- list(n = 1e4, R = 3, cvx = 0.01, cvy = 0.02, cil = 2, ciu = 10)
parameters_2 <- list(n = 1e4, R = 3, cvx = 1e-20, cvy = 0.02, cil = 2, ciu = 10, b0 = 0.03, b1 = 0.98)

# Structure tests
test_that(desc = "Check structure of simulate_eqa_data2()", code = {

  # Test data
  test_data_1 <- simulate_eqa_data2(parameters = parameters_1, type = 0L, AR = TRUE)

  # Check if list
  expect_true(object = is.list(test_data_1))

  # Check names
  expect_named(object = test_data_1,
               expected = c("SampleID", "ReplicateID", "MP_A", "MP_B"),
               ignore.order = FALSE,
               ignore.case = FALSE)

  # Check if lengths are equal
  expect_length(object = test_data_1$SampleID, n = 1e4 * 3)
  expect_length(object = test_data_1$ReplicateID, n = 1e4 * 3)
  expect_length(object = test_data_1$MP_A, n = 1e4 * 3)
  expect_length(object = test_data_1$MP_B, n = 1e4 * 3)

  # Check types
  expect_true(object = is.integer(test_data_1$SampleID))
  expect_true(object = is.integer(test_data_1$ReplicateID))
  expect_true(object = is.numeric(test_data_1$MP_A))
  expect_true(object = is.numeric(test_data_1$MP_B))

  # Check SampleID and ReplicateID
  expect_true(object = all(1:1e4 == unique(test_data_1$SampleID)))
  expect_true(object = all(1:3 == unique(test_data_1$ReplicateID)))
  expect_equal(object = test_data_1$SampleID[1:6], expected = rep(c(1, 2), each = 3))
  expect_equal(object = test_data_1$ReplicateID[1:6], expected = rep(1:3, 2))
})

# Value tests
test_that(desc = "Test if values correspond with given parameters", code = {

  # Test data
  test_data_1 <- simulate_eqa_data2(parameters = parameters_1, type = 0L, AR = TRUE)
  test_data_1_mor <- simulate_eqa_data2(parameters = parameters_1, type = 0L, AR = TRUE)
  test_data_2 <- simulate_eqa_data2(parameters = parameters_2, type = 0L, AR = FALSE)

  # Check cvx and cvy
  imprecision_estimates <- global_precision_estimates2(test_data_1)
  expect_true(object = abs(imprecision_estimates$Var_A - 0.0144) < 0.0017)
  expect_true(object = abs(imprecision_estimates$Var_B - 0.0036) < 0.0004)
  expect_true(object = abs(imprecision_estimates$CV_A - 0.02) < 0.0013)
  expect_true(object = abs(imprecision_estimates$CV_B - 0.01) < 0.0013)
  expect_true(object = abs(imprecision_estimates$lambda - 4) < 0.65)

  # Check cil and ciu
  expect_true(object = min(test_data_1_mor$MP_A) >= 2 - 5 * sqrt(0.0144 / 3))
  expect_true(object = min(test_data_1_mor$MP_B) >= 2 - 5 * sqrt(0.0036 / 3))
  expect_true(object = max(test_data_1_mor$MP_A) <= 10 + 5 * sqrt(0.0144 / 3))
  expect_true(object = max(test_data_1_mor$MP_B) <= 10 + 5 * sqrt(0.0036 / 3))

  # Check b0 and b1
  b1 <- var(test_data_2$MP_A, test_data_2$MP_B) / var(test_data_2$MP_B)
  b0 <- mean(test_data_2$MP_A) - b1 * mean(test_data_2$MP_B)
  expect_true(object = abs(b0 - 0.03) < 0.05)
  expect_true(object = abs(b1 - 0.98) < 0.005)

})



