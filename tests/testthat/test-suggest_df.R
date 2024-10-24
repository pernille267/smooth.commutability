library(fasteqa)
suppressWarnings(library(data.table))
library(testthat)

# Read test data
test_data <- fread(file = "~/Packages/test data smooth.commutability/fictive_crp_cs_data.csv")

# Testing ...
comparisons <- unique(test_data$comparison)
comp_id <- sample.int(10, size = 1)
complete_output <- suggest_df(data = test_data[comparison == comparisons[comp_id]], ns_score = TRUE, plots = TRUE, models = TRUE) |> suppressWarnings()
output_without_ns <- suggest_df(data = test_data[comparison == comparisons[comp_id]], ns_score = FALSE, plots = TRUE, models = TRUE)  |> suppressWarnings()
output_without_plots <- suggest_df(data = test_data[comparison == comparisons[comp_id]], ns_score = FALSE, plots = FALSE, models = TRUE)  |> suppressWarnings()
output_without_plots_and_models <- suggest_df(data = test_data[comparison == comparisons[comp_id]], ns_score = TRUE, plots = FALSE, models = FALSE)  |> suppressWarnings()
output_without_ns_plots_and_models <- suggest_df(data = test_data[comparison == comparisons[comp_id]], ns_score = FALSE, plots = FALSE, models = FALSE)  |> suppressWarnings()

test_that("suggest_df handles valid input correctly", {
  expect_s3_class(output_without_ns_plots_and_models, "data.table")
  expect_equal(nrow(output_without_ns_plots_and_models), 2)
  expect_equal(ncol(output_without_ns_plots_and_models), 7)
  expect_true(all(c("par", "zm_ar", "zd_ar", "zm_mor", "zd_mor", "loo", "gcv") %in% colnames(output_without_ns_plots_and_models)))
})

test_that("suggest_df handles ns_score parameter correctly", {
  expect_equal(nrow(output_without_plots$table), 2)
  expect_equal(nrow(output_without_plots_and_models), 3)
  expect_equal(nrow(output_without_ns_plots_and_models), 2)
  expect_true("non-linearity score" %in% output_without_plots_and_models$par)
})

test_that("suggest_df returns models when requested", {
  expect_type(complete_output, "list")
  expect_true("models" %in% names(complete_output))
  expect_type(complete_output$models, "list")
  expect_equal(length(complete_output$models), 6)
})

suggest_dfs <- test_data[, suppressWarnings(suggest_df(.SD, ns_score = TRUE)), by = comparison]

test_that("Check whether dfs are in valid ranges", {
  expect_true(all(suggest_dfs[par == "df"]$zm_ar >= 2 & suggest_dfs[par == "df"]$zm_ar <= 25))
  expect_true(all(suggest_dfs[par == "df"]$zd_ar >= 2 & suggest_dfs[par == "df"]$zd_ar <= 25))
  expect_true(all(suggest_dfs[par == "df"]$zm_mor >= 2 & suggest_dfs[par == "df"]$zm_mor <= 15))
  expect_true(all(suggest_dfs[par == "df"]$zd_mor >= 2 & suggest_dfs[par == "df"]$zd_mor <= 15))
})

test_that("Check whether non-linearity scores are in valid ranges", {
  expect_true(all(suggest_dfs[par == "non-linearity score"]$zm_ar >= -0.25))
  expect_true(all(suggest_dfs[par == "non-linearity score"]$zm_mor >= -0.5))
  expect_true(all(suggest_dfs[par == "non-linearity score"]$zd_ar >= -0.25))
  expect_true(all(suggest_dfs[par == "non-linearity score"]$zd_mor >= -0.5))
})

test_that("Check whether zetas are in valid ranges", {
  expect_true(all(suggest_dfs[par == "zeta"]$zm_ar >= 0.67))
  expect_true(all(suggest_dfs[par == "zeta"]$zm_mor >= 0.50))
  expect_true(all(suggest_dfs[par == "zeta"]$zd_ar >= 0.67))
  expect_true(all(suggest_dfs[par == "zeta"]$zd_mor >= 0.50))
})


set.seed(99)

pars <- list(n = 25, R = 3, cil = 2, ciu = 11, cvx = 0.01, cvy = 0.01, cve = 0)
sim_dat <- simulate_eqa_data2(pars, type = 2, TRUE) |> setDT()
output <- suggest_df(data = sim_dat, ns_score = TRUE, plots = TRUE)
dfminmor <- output$plots$second_deriv$df_mor[which.max(output$plots$second_deriv$second_deriv_zeta_mor)]
dfminar <- output$plots$second_deriv$df_ar[which.max(output$plots$second_deriv$second_deriv_zeta_ar)]
dfmaxmor <- obtain_df_max(output$plots$second_deriv$df_mor, output$plots$second_deriv$second_deriv_zeta_mor)
dfmaxar <- obtain_df_max(output$plots$second_deriv$df_ar, output$plots$second_deriv$second_deriv_zeta_ar)

test_that("Check whether df is correctly constrained", {
  expect_true(output$table$zm_mor[1] > dfminmor & output$table$zm_mor[1] < 2.01 * dfmaxmor)
  expect_true(output$table$zm_ar[1] > dfminar & output$table$zm_ar[1] < 2.01 * dfmaxar)
  expect_true(output$table$zd_mor[1] > dfminmor & output$table$zd_mor[1] < 2.01 * dfmaxmor)
  expect_true(output$table$zd_ar[1] > dfminar & output$table$zd_ar[1] < 2.01 * dfmaxar)

})
