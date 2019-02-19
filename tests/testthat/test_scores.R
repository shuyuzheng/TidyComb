context("calculate_score")
library("TidyComb")
library("dplyr")
response <- read.csv("response_test.csv", stringsAsFactors = FALSE)
response.mat <- reshape2::acast(response, conc_r ~ conc_c,
                                value.var = "response")
response.mat.noise <- AddNoise(response.mat, method = "random")
zip <- reshape2::acast(response, conc_r ~ conc_c, value.var = "synergy_zip")

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

# FitDoseResponse
fit <- FitDoseResponse(drug.row)
test_that("FitDoseResponse should return a 'drc' object.", {
  expect_equal(class(fit), "drc")
})
test_that("CalculateZIP should reture correct ZIP synergy score", {
  expect_equal(CalculateZIP(response.mat = response.mat), zip)
})