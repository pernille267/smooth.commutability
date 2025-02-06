library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

# Repr.
set.seed(99)

# Test data
test_parameters <- list(n = 25, R = 3, cvx = 0.01, cvy = 0.01, cil = 2, ciu = 10)
test_data <- simulate_eqa_data2(parameters = test_parameters, type = 2, AR = FALSE)
test_ss <- smoothing_spline(test_data, df = 4)

# Testing calculate_r() method
test_that(desc = "Check calculate_r()", code = {

  # get relevant components
  BTB_diag <- diag(test_ss$BTB)
  Omega_diag <- diag(test_ss$penalty_matrix)
  BTB_diag_restricted <- BTB_diag[-c(1:2, (length(BTB_diag) - 2):length(BTB_diag))]
  Omega_diag_restricted <- Omega_diag[-c(1:2, (length(Omega_diag) - 2):length(Omega_diag))]

  # Actual output and expected output
  actual <- calculate_r(test_ss$BTB, test_ss$penalty_matrix)
  expected <- sum(BTB_diag_restricted) / sum(Omega_diag_restricted)

  # Test if equal
  expect_equal(object = actual, expected = expected)

})

# Testing calculate_Q() method
test_that(desc = "Check calculate_Q()", code = {

  # Get relevant components
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda

  # Actual output and expected output
  actual <- calculate_Q(lambda, BTB, Omega)
  expected <- solve(BTB + lambda * Omega)

  # Test if equal
  expect_equal(object = actual, expected = expected)

})

# Testing calculate_S() method
test_that(desc = "Check calculate_S()", code = {

  # Get relevant components
  B <- test_ss$B
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda
  weightz <- test_ss$weights

  # Actual output and expected output
  actual <- calculate_S(lambda, weightz, B, BTB, Omega)
  expected <- B %*% solve(BTB + lambda * Omega) %*% t(B)

  # Test if equal
  expect_equal(object = actual, expected = expected)

})

# Testing calculate_S2() method
test_that(desc = "Check calculate_S2()", code = {

  # Get relevant components
  B <- test_ss$B
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda
  weightz <- test_ss$weights
  Q_expected <- solve(BTB + lambda * Omega)

  # Actual output and expected output
  actual_1 <- calculate_S2(weightz, B, calculate_Q(lambda, BTB, Omega))
  actual_2 <- calculate_S2(weightz, B, Q_expected)
  expected <- B %*% Q_expected %*% t(B)

  # Test if equal
  expect_equal(object = actual_1, expected = expected)
  expect_equal(object = actual_2, expected = expected)

})

# Testing calculate_cov_beta() method
test_that(desc = "Check calculate_cov_beta()", code = {

  # Get relevant components
  B <- test_ss$B
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda
  weightz <- test_ss$weights
  Q_expected <- solve(BTB + lambda * Omega)

  # Actual output and expected output
  actual_1 <- calculate_cov_beta(BTB, calculate_Q(lambda, BTB, Omega))
  actual_2 <- calculate_cov_beta(BTB, Q_expected)
  expected <- Q_expected %*% BTB %*% Q_expected

  # Test if equal
  expect_equal(object = actual_1, expected = expected)
  expect_equal(object = actual_2, expected = expected)

})

# Testing calculate_pred_var() method
test_that(desc = "Check calculate_pred_var()", code = {

  # Get relevant components
  B <- test_ss$B
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda
  weightz <- test_ss$weights
  Q_expected <- solve(BTB + lambda * Omega)
  cov_beta_expected <- Q_expected %*% BTB %*% Q_expected

  # Actual output and expected output
  actual_1 <- calculate_pred_var(B, calculate_cov_beta(BTB, calculate_Q(lambda, BTB, Omega)))[,]
  actual_2 <- calculate_pred_var(B, calculate_cov_beta(BTB, Q_expected))[,]
  actual_3 <- calculate_pred_var(B, cov_beta_expected)[,]
  expected <- diag(B %*% cov_beta_expected %*% t(B))

  # Test if equal
  expect_equal(object = actual_1, expected = expected)
  expect_equal(object = actual_2, expected = expected)
  expect_equal(object = actual_3, expected = expected)

})

# Testing calculate_df() method
test_that(desc = "Check calculate_df()", code = {

  # Get relevant components
  B <- test_ss$B
  BTB <- test_ss$BTB
  Omega <- test_ss$penalty_matrix
  lambda <- test_ss$lambda
  weightz <- test_ss$weights

  # Actual output and expected output
  actual <- calculate_df(lambda, weightz, B, BTB, Omega)
  expected <- sum(diag(B %*% solve(BTB + lambda * Omega) %*% t(B)))

  # Test if equal
  expect_equal(object = actual, expected = expected)

})
