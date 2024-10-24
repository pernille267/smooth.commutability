library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

is_lower_4_banded <- function(matrix) {
  n <- nrow(matrix)
  for (i in 1:n) {
    for (j in 1:n) {
      if (j > i + 4 && matrix[i, j] != 0) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

test_data <- simulate_eqa_data2(parameters = list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, cve = 0), type = 2)

x_orig <- test_data$MP_B[order(test_data$MP_B)]
x_unit <- (x_orig - min(x_orig)) / diff(range(x_orig))
k_inte <- calculate_interior_knots(x_unit)
Omega <- penalty_matrix(x_min = 0, x_max = 1, interior_knots = k_inte)

test_that("Check properties of penalty matrix", {
  expect_true(isSymmetric(Omega), label = "[Penalty matrix is symmetric]")
  expect_true(is_lower_4_banded(Omega), label = "[Penalty matrix is lower 4-banded]")
  expect_true(all(abs(rowSums(Omega) - 0) < 1e-3), label = "[All row sums of penalty matrix = 0]")
  expect_true(all(diag(Omega) >= 0), label = "[All diagonal elements of penalty matrix are positive]")
  expect_true(all(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values >= -1e-3), label = "[Penalty matrix is non-negative definite]")
})


