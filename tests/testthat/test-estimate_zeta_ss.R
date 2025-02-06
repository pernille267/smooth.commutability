library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

# Initialize
detailed_testing <- TRUE

# Testing on simulated data

# Parameters used for testing
parameters_1 <- list(n = 40, R = 4, cil = 2, ciu = 10, cvx = 0.01, cvy = 0.01, b0 = 0.025, b1 = 0.975)
parameters_2 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 0.05, cvy = 0.01)
parameters_3 <- list(n = 25, R = 3, cil = 2, ciu = 10, cvx = 1e-20, cvy = 0.01, qpos = 1, qran = 0.5, mmax = 5)

# Tests 1.1.1 - 1.3.2
test_that("Check structure of output of estimate_zeta_ss()", code = {

  # Test data
  test_data_1 <- simulate_eqa_data2(parameters = parameters_1,
                                    type = 0,
                                    AR = TRUE)

  # Actual output
  actual_1 <- estimate_zeta_ss(data = test_data_1, df = 2)
  actual_2 <- estimate_zeta_ss(data = test_data_1, df = 5)

  ## Check names ---------------------------------------------------------------
  expect_named(actual_1, "zeta",
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.1.1]")
  expect_named(actual_2, "zeta",
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 1.1.2]")

  ## Check class ---------------------------------------------------------------
  expect_true(is.list(actual_1),
              label = "[TEST 1.2.1]")
  expect_true(is.list(actual_2),
              label = "[TEST 1.2.2]")

  ## Check type ----------------------------------------------------------------
  expect_true(is.double(actual_1$zeta),
              label = "[TEST 1.3.1]")
  expect_true(is.double(actual_2$zeta),
              label = "[TEST 1.3.2]")
})

# Tests 2.1.1 - 2.2.5
test_that("Check value of estimate_zeta_ss()", code = {

  # Test data
  test_data_1 <- simulate_eqa_data2(parameters = parameters_1,
                                    type = 0,
                                    AR = TRUE)
  test_data_2 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 1,
                                    AR = TRUE)
  test_data_3 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 2,
                                    AR = TRUE)
  test_data_4 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 3,
                                    AR = TRUE)
  test_data_5 <- simulate_eqa_data2(parameters = parameters_3,
                                    type = 2,
                                    AR = TRUE)

  # Actual output (I)
  actual_1 <- estimate_zeta_ss(data = test_data_1)
  actual_2 <- estimate_zeta_ss(data = test_data_2)
  actual_3 <- estimate_zeta_ss(data = test_data_3)
  actual_4 <- estimate_zeta_ss(data = test_data_4)
  actual_5 <- estimate_zeta_ss(data = test_data_5)

  # Check if estimated zeta >= 0.5
  expect_gte(actual_1$zeta,
             expected = 0.5,
             label = "[TEST 2.1.1]")
  expect_gte(actual_2$zeta,
             expected = 0.5,
             label = "[TEST 2.1.2]")
  expect_gte(actual_3$zeta,
             expected = 0.5,
             label = "[TEST 2.1.3]")
  expect_gte(actual_4$zeta,
             expected = 0.5,
             label = "[TEST 2.1.4]")
  expect_gte(actual_5$zeta,
             expected = 0.5,
             label = "[TEST 2.1.5]")

  # Actual output (II)
  actual_6 <- estimate_zeta_ss(data = test_data_1, mor = TRUE)
  actual_7 <- estimate_zeta_ss(data = test_data_2, mor = TRUE)
  actual_8 <- estimate_zeta_ss(data = test_data_3, mor = TRUE)
  actual_9 <- estimate_zeta_ss(data = test_data_4, mor = TRUE)
  actual_0 <- estimate_zeta_ss(data = test_data_5, mor = TRUE)

  # Check if estimated zeta >= 0
  expect_gte(actual_1$zeta,
             expected = 0,
             label = "[TEST 2.2.1]")
  expect_gte(actual_2$zeta,
             expected = 0,
             label = "[TEST 2.2.2]")
  expect_gte(actual_3$zeta,
             expected = 0,
             label = "[TEST 2.2.3]")
  expect_gte(actual_4$zeta,
             expected = 0,
             label = "[TEST 2.2.4]")
  expect_gte(actual_5$zeta,
             expected = 0,
             label = "[TEST 2.2.5]")

})

# Tests 3.1.1 - 3.1.5
test_that("Check value of estimate_zeta_ss() when df <= 2", code = {

  # Test data
  test_data_1 <- simulate_eqa_data2(parameters = parameters_1,
                                    type = 0,
                                    AR = TRUE)
  test_data_2 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 1,
                                    AR = TRUE)
  test_data_3 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 2,
                                    AR = TRUE)
  test_data_4 <- simulate_eqa_data2(parameters = parameters_2,
                                    type = 3,
                                    AR = TRUE)
  test_data_5 <- simulate_eqa_data2(parameters = parameters_3,
                                    type = 2,
                                    AR = TRUE)

  # Actual output
  actual_1 <- estimate_zeta_ss(data = test_data_1, df = 2)
  actual_2 <- estimate_zeta_ss(data = test_data_2, df = 1)
  actual_3 <- estimate_zeta_ss(data = test_data_3, df = 0)
  actual_4 <- estimate_zeta_ss(data = test_data_4, df = 1)
  actual_5 <- estimate_zeta_ss(data = test_data_5, df = 2)

  expect_equal(object = actual_1$zeta,
               expected = estimate_zeta(test_data_1)$zeta,
               tolerance = 1e-4,
               label = "[TEST 3.1.1]")
  expect_equal(object = actual_2$zeta,
               expected = estimate_zeta(test_data_2)$zeta,
               tolerance = 1e-4,
               label = "[TEST 3.1.2]")
  expect_equal(object = actual_3$zeta,
               expected = estimate_zeta(test_data_3)$zeta,
               tolerance = 1e-4,
               label = "[TEST 3.1.3]")
  expect_equal(object = actual_4$zeta,
               expected = estimate_zeta(test_data_4)$zeta,
               tolerance = 1e-4,
               label = "[TEST 3.1.4]")
  expect_equal(object = actual_5$zeta,
               expected = estimate_zeta(test_data_5)$zeta,
               tolerance = 1e-4,
               label = "[TEST 3.1.5]")
})


# Testing on "real data"
real_data <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

# Tests 4.1.1 - 4.3.3
test_that(desc = "Check value of estimate_zeta_ss() (II)", code = {

  # Test data
  real_data_1 <- real_data[comparison == "AQT90 - Chroma"]
  real_data_2 <- real_data[comparison == "Luminex - QuickRead"]
  real_data_3 <- real_data[comparison == "Chroma - QuickRead"]

  # Actual output
  actual_1 <- estimate_zeta_ss(data = real_data_1, df = 2)
  actual_2 <- estimate_zeta_ss(data = real_data_1, df = 3)
  actual_3 <- estimate_zeta_ss(data = real_data_1, df = 6)
  actual_4 <- estimate_zeta_ss(data = real_data_2, df = 2)
  actual_5 <- estimate_zeta_ss(data = real_data_2, df = 3)
  actual_6 <- estimate_zeta_ss(data = real_data_2, df = 6)
  actual_7 <- estimate_zeta_ss(data = real_data_3, df = 2)
  actual_8 <- estimate_zeta_ss(data = real_data_3, df = 3)
  actual_9 <- estimate_zeta_ss(data = real_data_3, df = 6)

  # Check whether class = list
  expect_true(object = is.list(actual_2),
              label = "[TEST 4.1.1]")
  expect_true(object = is.list(actual_5),
              label = "[TEST 4.1.2]")
  expect_true(object = is.list(actual_9),
              label = "[TEST 4.1.3]")

  # Check whether zeta is numeric and strictly positive
  expect_true(object = is.numeric(actual_3$zeta) & actual_3$zeta >= 0.0 & actual_3$zeta <= 1e3,
              label = "[TEST 4.2.1]")
  expect_true(object = is.numeric(actual_6$zeta) & actual_6$zeta >= 0.0 & actual_2$zeta <= 1e3,
              label = "[TEST 4.2.2]")
  expect_true(object = is.numeric(actual_8$zeta) & actual_8$zeta >= 0.0 & actual_8$zeta <= 1e3,
              label = "[TEST 4.2.3]")

  # Check whether zeta values are plausible for df = 2
  impr_1 <- global_precision_estimates2(real_data_1)
  impr_2 <- global_precision_estimates2(real_data_2)
  impr_3 <- global_precision_estimates2(real_data_3)
  lambda_1 <- impr_1$lambda
  lambda_2 <- impr_2$lambda
  lambda_3 <- impr_3$lambda
  sigma_h_squared_1 <- sigma_h_squared_deming2(data = real_data_1, lambda = lambda_1)
  sigma_h_squared_2 <- sigma_h_squared_deming2(data = real_data_2, lambda = lambda_2)
  sigma_h_squared_3 <- sigma_h_squared_deming2(data = real_data_3, lambda = lambda_3)
  zeta_1 <- (sigma_h_squared_1 + lambda_1 * sigma_h_squared_1) / (impr_1$Var_A + impr_1$Var_B)
  zeta_2 <- (sigma_h_squared_2 + lambda_2 * sigma_h_squared_2) / (impr_2$Var_A + impr_2$Var_B)
  zeta_3 <- (sigma_h_squared_3 + lambda_3 * sigma_h_squared_3) / (impr_3$Var_A + impr_3$Var_B)
  expect_true(object = actual_1$zeta >= zeta_1 * 0.75 & actual_1$zeta <= zeta_1 * 1.25,
              label = "[TEST 4.3.1]")
  expect_true(object = actual_4$zeta >= zeta_2 * 0.75 & actual_4$zeta <= zeta_2 * 1.25,
              label = "[TEST 4.3.2]")
  expect_true(object = actual_7$zeta >= zeta_3 * 0.75 & actual_7$zeta <= zeta_3 * 1.25,
              label = "[TEST 4.3.3]")

})

# Tests 5.1.1 - 5.4.3
test_that(desc = "Check structure of estimate_zetas_ss()", code = {

  # Actual output
  actual_1 <- estimate_zetas_ss(data = real_data,
                                df = NULL,
                                weighted = FALSE)
  actual_2 <- estimate_zetas_ss(data = real_data,
                                df = 5,
                                weighted = FALSE)
  actual_3 <- estimate_zetas_ss(data = real_data,
                                df = runif(10, 2, 10),
                                weighted = FALSE)

  # Check if data.table class
  expect_true(object = is.data.table(actual_1),
              label = "[TEST 5.1.1]")
  expect_true(object = is.data.table(actual_2),
              label = "[TEST 5.1.2]")
  expect_true(object = is.data.table(actual_3),
              label = "[TEST 5.1.3]")

  # Check number of rows and number of columns
  expect_true(object = nrow(actual_1) == 10L && ncol(actual_1) == 2L,
              label = "[TEST 5.2.1]")
  expect_true(object = nrow(actual_2) == 10L && ncol(actual_2) == 2L,
              label = "[TEST 5.2.2]")
  expect_true(object = nrow(actual_3) == 10L && ncol(actual_3) == 2L,
              label = "[TEST 5.2.3]")

  # Check column names
  expect_named(object = actual_1,
               expected = c("comparison", "zeta"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 5.3.1]")
  expect_named(object = actual_2,
               expected = c("comparison", "zeta"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 5.3.2]")
  expect_named(object = actual_3,
               expected = c("comparison", "zeta"),
               ignore.order = FALSE,
               ignore.case = FALSE,
               label = "[TEST 5.3.3]")

  # Check comparison column
  expected_comparisons <- unique(real_data$comparison)
  expect_true(object = all(actual_1$comparison == expected_comparisons),
              label = "[TEST 5.4.1]")
  expect_true(object = all(actual_2$comparison == expected_comparisons),
              label = "[TEST 5.4.2]")
  expect_true(object = all(actual_3$comparison == expected_comparisons),
              label = "[TEST 5.4.3]")

})

# Tests 6.1.1 - 6.3.2
test_that(desc = "Check values of estimate_zetas_ss()", code = {

  # Test data
  real_data_1 <- real_data[comparison == "AQT90 - Chroma"]
  real_data_2 <- real_data[comparison == "Luminex - QuickRead"]
  real_data_3 <- real_data[comparison == "Chroma - QuickRead"]

  # Actual output
  random_dfs <- runif(10, 2, 10)
  actual_1 <- estimate_zetas_ss(data = real_data,
                                df = NULL,
                                weighted = FALSE)
  actual_2 <- estimate_zetas_ss(data = real_data,
                                df = 5,
                                weighted = FALSE)
  actual_3 <- estimate_zetas_ss(data = real_data,
                                df = random_dfs,
                                weighted = FALSE)

  actual_4 <- estimate_zeta_ss(data = real_data_1)
  actual_5 <- estimate_zeta_ss(data = real_data_2, df = 5)
  actual_6 <- estimate_zeta_ss(data = real_data_3, df = random_dfs[7])

  # Check if zeta is numeric
  expect_true(object = is.numeric(actual_1$zeta),
              label = "[TEST 6.1.1]")
  expect_true(object = is.numeric(actual_2$zeta),
              label = "[TEST 6.1.2]")
  expect_true(object = is.numeric(actual_3$zeta),
              label = "[TEST 6.1.3]")

  # Check if all zeta values are >= 0
  expect_true(object = all(actual_1$zeta >= 0.0),
              label = "[TEST 6.2.1]")
  expect_true(object = all(actual_2$zeta >= 0.0),
              label = "[TEST 6.2.2]")
  expect_true(object = all(actual_3$zeta >= 0.0),
              label = "[TEST 6.2.3]")

  # Check if zeta values from estimate_zetas_ss() correspond with estimate_zeta_ss()
  expect_true(object = abs(actual_1[comparison == "AQT90 - Chroma"]$zeta - actual_4$zeta) <= 1e-6,
              label = "[TEST 6.3.1]")
  expect_true(object = abs(actual_2[comparison == "Luminex - QuickRead"]$zeta - actual_5$zeta) <= 1e-6,
              label = "[TEST 6.3.2]")
  expect_true(object = abs(actual_3[comparison == "Chroma - QuickRead"]$zeta - actual_6$zeta) <= 1e-6,
              label = "[TEST 6.3.3]")

})

# Tests 7.1.1 - 7.2.6
test_that(desc = "Error handling", code = {

  # Test data
  real_data_1 <- real_data[comparison == "AQT90 - Chroma"]
  real_data_2 <- real_data[comparison == "Luminex - QuickRead"]
  real_data_3 <- real_data[comparison == "Chroma - QuickRead"]
  names(real_data_3) <- c("comparison", "SID", "RID", "MP_A", "MP_B")

  # Error handling for estimate_zeta_ss()
  expect_error(object = estimate_zeta_ss(data = real_data_1, df = 100),
               regexp = "df cannot exceed the number of unique samples",
               label = "[TEST 7.1.1]")
  expect_error(object = estimate_zeta_ss(data = real_data_1, df = "wrong!"),
               regexp = "but is a character",
               label = "[TEST 7.1.2]")
  expect_error(object = estimate_zeta_ss(data = real_data_2, df = FALSE),
               regexp = "but is a logical",
               label = "[TEST 7.1.3]")
  expect_error(object = estimate_zeta_ss(data = real_data_1, weighted = 3),
               regexp = "but is a numeric",
               label = "[TEST 7.1.4]")
  expect_error(object = estimate_zeta_ss(data = real_data_2, weighted = NULL),
               regexp = "but is NULL",
               label = "[TEST 7.1.5]")
  expect_error(object = estimate_zeta_ss(data = real_data_3),
               regexp = "does not have the necessary columns",
               label = "[TEST 7.1.6]")

  # Error handling for estimate_zetas_ss()
  # ...

})


# These tests takes a long time, so avoid it...
# These tests are also approximate, so they will only trigger if something
# is very wrong!
if(detailed_testing){

  test_that(desc = "Testing relationship between MOR and AR", code = {
    mor_and_ar <- sapply(X = 1:1000,
                         FUN = function(x){
                          test_data_2 <- simulate_eqa_data2(parameters = parameters_2,
                                                            type = 2,
                                                            AR = TRUE)
                          return(list("AR" = estimate_zeta_ss(data = test_data_2, mor = FALSE)$zeta,
                                      "MOR" = estimate_zeta_ss(data = test_data_2, mor = TRUE)$zeta))
                        }, simplify = FALSE)
    mor_and_ar <- rbindlist(mor_and_ar)

    # Check correlation
    expect_gte(object = cor(mor_and_ar$AR, mor_and_ar$MOR),
               expected = 0.50,
               label = "[TEST D.1.1]")

    # Check ratio of standard deviations
    expect_gte(object = sd(mor_and_ar$MOR, na.rm = TRUE) / sd(mor_and_ar$AR, na.rm = TRUE),
               expected = 2.43,
               label = "[TEST D.1.2]")

    # Check median of ratio
    median_mor_ar_ratio <- median(mor_and_ar$MOR / mor_and_ar$AR, na.rm = TRUE)
    expect_true(object = median_mor_ar_ratio >= 1.00 & median_mor_ar_ratio <= 1.25,
                label = "[TEST D.1.3]")

    # Check mad of ratio
    mad_mor_ar_ratio <- mad(mor_and_ar$MOR / mor_and_ar$AR, na.rm = TRUE)
    lower_limit <- 1/3 - 1/ 7 / 3
    upper_limit <- 1/3 + 1/ 7 / 3
    expect_true(object = mad_mor_ar_ratio >= lower_limit & mad_mor_ar_ratio <= upper_limit,
                label = "[TEST D.1.4]")

  })
}







