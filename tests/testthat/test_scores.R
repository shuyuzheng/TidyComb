context("calculate_score")
library("TidyComb")
library("dplyr")
response <- read.csv("response_with_scores.csv", stringsAsFactors = FALSE)
response.mat <- reshape2::acast(response, conc_r ~ conc_c,
                                value.var = "response")

# ExtractSingleDrug
drug.row <- ExtractSingleDrug(response.mat = response.mat, "row")
drug.col <- ExtractSingleDrug(response.mat = response.mat, "col")
test_that("ExtractSingleDrug should return a matrix with two columns 'conc' and
          'effect'.", {
  expect_equal(colnames(drug.row), c("response", "dose"))
  expect_equal(colnames(drug.col), c("response", "dose"))
})

test_that("ExtractSingleDrug extract correct information of drug on correct
          dim", {
  expect_equal(drug.row, response %>% filter(conc_c == 0) %>%
                           select(response, dose = conc_r) %>% arrange(dose))
  expect_equal(drug.col, response %>% filter(conc_r == 0) %>%
                           select(response, dose = conc_c) %>% arrange(dose))
})

test_that("Value passed to ExtractSingleDrug dim should be either 'row' or
          'col'", {
  expect_error(ExtractSingleDrug(response.mat = response.mat, "whatever"),
               "Values for 'dim' should be eighther 'row' or 'col'!")
})

# FittingSingleDrug
fit <- FittingSingleDrug(drug.row)
test_that("FittingSingleDrug should return a list with 2 element 'fitted' and
          'model'", {
  expect_equal(class(fit), "list")
  expect_equal(names(fit), c("fitted", "model"))
  expect_equal(class(fit$model), "drc")
  expect_equal(class(fit), "numeric")
})
test_that("CalculateZIP should reture correct ZIP synergy score", {
  expect_equal(CalculateZIP(response.mat = response.mat), response$synergy_zip)
})