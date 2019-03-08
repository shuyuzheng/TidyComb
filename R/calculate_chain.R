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
#'       \item Calculate Synergy Scores
#'       \item Calculate CSS
#'     }
#'   \item Summarize and generate surface
#' }
#' @param response.mat
#' @param correction
#' @param ...
#'
#' @return
#' @export
CalculateMat <- function(response.mat, correction = TRUE, ...) {
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
  col.type <- FindModelType(col.type)

  # drug_row
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  row.model <- FitDoseResponse(drug.row)
  row.type <- FindModelType(row.type)

  # 2.2 Extract coeficients and IC50
  # drug_col
  col.coe <- stats::coef(col.model)
  # drug_rowr
  row.coe <- stats::coef(row.model)

  # 2.3 Calculate DSS (using single drug response but without that at 0
  # concentration)
  col.dss <- CalculateSens(drug.col)
  row.dss <- CalculateSens(drug.row)

  # 3. whole matrix process
  # 3.1 Calculate synergyscores
  zip <- CalculateZIP(response.mat, drug.row.model = row.model,
                      drug.col.model = col.model)
  loewe <- CalculateLoewe(response.mat, drug.row.type = row.type,
                          drug.row.par = row.coe, drug.col.type = col.type,
                          drug.col.par = col.coe)
  hsa <- CalculateHSA(response.mat)
  bliss <- CalculateBliss(response.mat)

  synergy <- Reduce(function(x, y) merge(x = x, y = y,
                                         by = c("conc_c", "conc_r")),
                    list(zip, loewe, hsa, bliss))
  # 3.2 Calculate CSS
  col.ic50 <- CalculateIC50(col.coe, col.type, max(drug.col$dose))
  row.ic50 <- CalculateIC50(row.coe, row.type, max(drug.row$dose))

  res <- ImputeIC50(response.mat, col.ic50 = col.ic50, row.ic50 = row.ic50)
  # a particular row selected according to ic50_row
  col.css <- CalculateSens(res$tempcf_c)
  # a particular column selected according to ic50_col
  row.css <- CalculateSens(res$tempcf_r)

  css <- mean(col.css, row.css)

  S <- css - mean(col.dss, row.dss)
}