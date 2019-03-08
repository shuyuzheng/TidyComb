# TidyComb
# Functions for pipe all calculation togather
# Copyright Shuyu Zheng
#
# Functions on this page:
#

#' Chain all calculatrion about one dose-response matrix together
#'
#' \code{CalculateMat} chains all calculations about one dose-response matrix (
#' one drug-drug interaction block) together. The calculations includes:
#' dose-response curve fitting, synergy scores (ZIP, Bliss, Loewe,
#' HSA), drug sensitivity (DSS, CSS), generate drug-drug response surface and
#' generating summary scores for each block.
#'
#' The steps for calculation:
#' \enumerate{
#'   \item Pre-process Matrix
#'     \enumerate{
#'       \item Impute for missing value in original matrix.
#'       \item Add noise to original matrix.
#'       \item Correct baseline to 0, if \code{correction} is \code{TRUE}.
#'     }
#'   \item Single drug process
#'     \enumerate{
#'       \item Extract and fitting single drugs.
#'       \item Extract coeficients from fitted model. (b, c, d, e, IC50)
#'       \item Calculate DSS
#'     }
#'   \item Whole response matrix process
#'     \enumerate{
#'       \item Fit dose-response model through whole matrix
#'       \item Calculate CSS
#'       \item Calculate Synergy Scores
#'     }
#'   \item Summarize and generate surface
#' }
#' @param response.mat
#' @param correction
#' @param ...
#'
#' @return
#' @export
CalculateMat <- function(response.mat, correction = TRUE, ...){
  # 1. Pre-processing
  # 1.1 Impute for missing value in original matrix
  response.mat <- ImputeNear(response.mat, times = 2)

  # 1.2. Add random noise to original matrix
  response.mat <- AddNoise(response.mat, method = "random")

  # 1.3. Correct baseline to 0 if "correction" is TRUE
  if (correction) {
    response.mat <- CorrectBaseLine(response.mat)
  }

  # 2. Single drug process
  # 2.1. Fit single drug dose-response curve
  # drug_col
  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  col.model <- FitDoseResponse(drug.col)
  col.method <- FindModelType(col.model)

  # drug_row
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  row.model <- FitDoseResponse(drug.row)
  row.method <- FindModelType(row.model)

  # 2.2 Extract coeficients and IC50
  # drug_col
  col.coe <- ExtractCoefs(col.model, col.method, max(drug.col$dose))
  # drug_rowr
  row.coe <- ExtractCoefs(row.model, row.method, max(drug.row$dose))

  # 2.3 Calculate DSS


}