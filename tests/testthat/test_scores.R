context("calculate_score")
library("TidyComb")
library("dplyr")
load("synergy.rda")

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
                 dplyr::select(response, dose = conc_r) %>% arrange(dose))
  expect_equal(drug.col, response %>% filter(conc_r == 0) %>%
                 dplyr::select(response, dose = conc_c) %>% arrange(dose))
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

# BaselineCorrect
correct.response <- CorrectBaseLine(response.mat.noise)
# Synergy scores
test_that("CalculateZIP should reture correct ZIP synergy score", {
  expect_equal(CalculateZIP(response.mat = correct.response), zip)
})

test_that("CalculateHSA should reture correct HSA synergy score", {
  expect_equal(CalculateHSA(response.mat = correct.response), hsa)
})

test_that("CalculateBliss should reture correct Bliss synergy score", {
  expect_equal(CalculateBliss(response.mat = correct.response), bliss)
})

test_that("CalculateLoewe should reture correct Loewe synergy score", {
  expect_equal(CalculateLoewe(response.mat = correct.response), loewe)
})