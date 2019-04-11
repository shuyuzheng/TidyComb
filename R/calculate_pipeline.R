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
#'       \item Calculate Surface
#'       \item Calculate CSS
#'     }
#'   \item Summarize and generate surface
#' }
#' @param response.mat A matrix contain the drug combination reaponse value.
#' Its column names are doses of drug added along columns. Its row name are
#' doses of drug added along rows.
#' @param correction A logic value to indicate whether or not do baseline
#' correction before calculation
#' @param ... Other argumants required by nested functions
#'
#' @return It contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified response value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, dss, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed response value and
#'     synergy scores of each drug dose response pair.
#'  }
#' @export
CalculateMat <- function(response.mat, correction = TRUE, ...) {

  options(scipen = 999)

  response.out <- response.mat

  # 1. Pre-processing
  # 1.1 Impute for missing value in original matrix
  response.mat <- ImputeNear(response.mat, times = 2)

  # 1.3. Correct baseline to 0 if "correction" is TRUE
  if (correction) {
    response.mat <- CorrectBaseLine(response.mat)
  }
  # 2. Single drug process
  # 2.1. Fit single drug dose-response curve
  # drug_col
  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  col.model <- FitDoseResponse(drug.col)
  col.type <- as.character(FindModelType(col.model))

  # drug_row
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  row.model <- FitDoseResponse(drug.row)
  row.type <- as.character(FindModelType(row.model))

  # 2.2 Extract coeficients
  # drug_col
  col.coe <- stats::coef(col.model)
  # drug_rowr
  row.coe <- stats::coef(row.model)

  curves <- as.data.frame(rbind(col.coe, row.coe))
  colnames(curves) <- sub(":(Intercept)", "", colnames(curves), fixed = TRUE)
    curves$model <- c(col.type, row.type)
  curves$dim <- c("col", "row")

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

  synergy <- lapply(list(response.out, zip, loewe, hsa, bliss), reshape2::melt)
  synergy <- Reduce(function(x, y) {
    merge(x = x, y = y, by = c("Var1", "Var2"))}, synergy)

  colnames(synergy) <- c("conc_r", "conc_c", "response", "synergy_zip",
                         "synergy_loewe", "synergy_hsa", "synergy_bliss")

  # 3.2 calculate surface
  smooth.res <- smoothing(response.mat)
  smooth.zip <- smoothing(zip)
  smooth.hsa <- smoothing(hsa)
  smooth.bliss <- smoothing(bliss)
  smooth.loewe <- smoothing(loewe)

  surface <- lapply(list(smooth.res, smooth.zip, smooth.loewe, smooth.hsa,
                         smooth.bliss), reshape2::melt)
  surface <- Reduce(function(x, y) {
    merge(x = x, y = y, by = c("Var1", "Var2"))}, surface)

  colnames(surface) <- c("conc_r", "conc_c", "response", "synergy_zip",
                         "synergy_loewe", "synergy_hsa", "synergy_bliss")

  # 3.3 Calculate CSS
  col.ic50 <- CalculateIC50(col.coe, col.type, max(drug.col$dose))
  row.ic50 <- CalculateIC50(row.coe, row.type, max(drug.row$dose))

  imputed.ic50 <- ImputeIC50(response.mat, col.ic50 = col.ic50,
                             row.ic50 = row.ic50)
  # a particular row selected according to ic50_row
  col.css <- CalculateSens(imputed.ic50$tempcf_c)
  # a particular column selected according to ic50_col
  row.css <- CalculateSens(imputed.ic50$tempcf_r)

  css <- mean(col.css, row.css)

  S <- css - max(col.dss, row.dss)

  sum <- synergy %>%
    dplyr::filter(conc_r != 0 & conc_c != 0) %>%
    dplyr::select(-conc_r, -conc_c) %>%
    apply(2, mean, rm.na = TRUE)

  sum <- data.frame(t(sum), ic50_row = row.ic50 , ic50_col = col.ic50,
                    dss_row = row.dss, dss_col = col.dss, css_row = row.css,
                    css_col = col.css, css = css, S = S)

  res <- list(synergy = synergy, surface = surface,
              summary = sum, curve = curves)
  return(res)

  # clean up
  gc()
}

#' Calculate Drug Combination data in template format
#'
#' @param template a dataframe in the format as template. Columns "block_id",
#' "drug_row", "drug_col", "response", "conc_r", "conc_c", "conc_r_unit",
#' "conc_c_unit","cell_line_name", "drug_row", "drug_col" are reqired.
#'
#' @param ... Other arguments required by nested functions
#'
#' @return A list. It contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified response value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, dss, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed response value and
#'     synergy scores of each drug dose response pair.
#'  }
#'
#' @export
CalculateTemplate <- function(template, ...) {
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
             "conc_r_unit", "conc_c_unit","cell_line_name", "drug_row",
             "drug_col") %in%
           colnames(template)))
    stop("The input data must contain the following columns: ",
         "block_id, drug_row, drug_col, response,\n",
         "conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
         "cell_line_name, drug_row, drug_col.")

  blocks <- unique(template$block_id)

  # generate container
  synergy <- data.frame(block_id = integer(), conc_r = numeric(),
                        conc_c = numeric(), response = numeric(),
                        synergy_zip = numeric(), synergy_bliss = numeric(),
                        synergy_loewe = numeric(), synergy_hsa = numeric(),
                        stringsAsFactors = FALSE)
  surface <- synergy
  curve <- data.frame(block_id = integer(), b = numeric(),
                      c = numeric(), d = numeric(), e = numeric(),
                      model = numeric(), drug.row = numeric(),
                      drug.col = numeric, stringsAsFactors = FALSE)
  summary <- data.frame(block_id = integer(), synergy_zip = numeric(),
                        synergy_bliss = numeric(), synergy_hsa = numeric(),
                        synergy_loewe = numeric(), ic50_row = numeric() ,
                        ic50_col = numeric(), dss_row = numeric(),
                        dss_col = numeric(), css_row = numeric(),
                        css_col = numeric(), css = numeric(), S = numeric(),
                        stringsAsFactors = FALSE)

  for (block in blocks) {
    # 1. Generate response matrix for each block
    response <- template %>%
      dplyr::filter(block_id == block)

    # 1.2. Add random noise to original matrix
    response <- AddNoise(response, method = "random")
    response.mat <- response %>%
      dplyr::select(conc_r, conc_c, response) %>%
      reshape2::acast(conc_r ~ conc_c, value.var = "response")

    # 2. Do calculation on matrix
    tmp <- CalculateMat(response.mat = response.mat)
    tmp <- lapply(tmp, function(x){
      x$block_id = rep(block, nrow(x))
      return(x)
      })

    # 3. Add information to summary table
    info <- response %>%
      dplyr::select(drug_row, drug_col, cell_line_name,
                    conc_r_unit, conc_c_unit) %>%
      unique()

    if (nrow(info) > 1) {
      warning("The summary data of block ", block, " contains ", nrow(info),
              " different versions.")
    } else if (nrow(info) < 1) {
      warning("The summary data of block ", block, " is missing.")
    }

    tmp$summary <- cbind.data.frame(info, tmp$summary)

    # 4. fill drug names to curve table

    tmp$curve$drug_col <- rep(NA, 2)
    tmp$curve$drug_col[which(tmp$curve$dim ==
                               "col")] <- as.character(info$drug_col)

    tmp$curve$drug_row <- rep(NA, 2)
    tmp$curve$drug_row[which(tmp$curve$dim ==
                               "row")] <- as.character(info$drug_row)

    # 4. Append new tables to containers
    synergy <- rbind.data.frame(synergy, tmp$synergy)
    surface <- rbind.data.frame(surface, tmp$surface)
    curve <- rbind.data.frame(curve, tmp$curve)
    summary <- rbind.data.frame(summary, tmp$summary)

    # Clean temporary file
    tmp <- list()
    info <- data.frame()
    response <- data.frame()
    respons.mat <- matrix()
  }

  curve <- dplyr::select(curve, block_id, drug_row, drug_col, b, c, d, e, model)

  return(list(synergy = synergy, surface = surface,
              curve = curve, summary = summary))
}

multiResultClass <- function(synergy=NULL, summary=NULL, surface = NULL,
                             curve = NULL) {
  me <- list(
    synergy = synergy,
    summary = summary,
    surface = surface,
    curve = curve
  )

  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}

#' Parallel Calculate Drug Combination data in template format
#'
#' @param template a dataframe in the format as template. Columns "block_id",
#' "drug_row", "drug_col", "response", "conc_r", "conc_c", "conc_r_unit",
#' "conc_c_unit","cell_line_name", "drug_row", "drug_col" are reqired.
#'
#' @param cores A integer. It indicates number of cores would be allocated to
#' the parallel processed
#'
#' @param ... Other arguments required by nested functions
#'
#' @return A list. It contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified response value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, dss, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed response value and
#'     synergy scores of each drug dose response pair.
#'  }
#'
#' @export
ParCalculateTemplate <- function(template, cores = 1, ...) {
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
             "conc_r_unit", "conc_c_unit","cell_line_name", "drug_row",
             "drug_col") %in%
           colnames(template)))
    stop("The input data must contain the following columns: ",
         "block_id, drug_row, drug_col, response,\n",
         "conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
         "cell_line_name, drug_row, drug_col.")

  blocks <- unique(template$block_id)

  # generate container
  # synergy <- data.frame(block_id = integer(), conc_r = numeric(),
  #                       conc_c = numeric(), response = numeric(),
  #                       synergy_zip = numeric(), synergy_bliss = numeric(),
  #                       synergy_loewe = numeric(), synergy_hsa = numeric(),
  #                       stringsAsFactors = FALSE)
  # surface <- synergy
  # curve <- data.frame(block_id = integer(), b = numeric(),
  #                     c = numeric(), d = numeric(), e = numeric(),
  #                     model = numeric(), drug.row = numeric(),
  #                     drug.col = numeric, stringsAsFactors = FALSE)
  # summary <- data.frame(block_id = integer(), synergy_zip = numeric(),
  #                       synergy_bliss = numeric(), synergy_hsa = numeric(),
  #                       synergy_loewe = numeric(), ic50_row = numeric() ,
  #                       ic50_col = numeric(), dss_row = numeric(),
  #                       dss_col = numeric(), css_row = numeric(),
  #                       css_col = numeric(), css = numeric(), S = numeric(),
  #                       stringsAsFactors = FALSE)

  registerDoParallel(cores)

  res <- foreach::foreach (i = 1:length(blocks)) %dopar% {
    # 1. Generate response matrix for each block
    result <- multiResultClass()
    response <- template %>%
      dplyr::filter(block_id == blocks[i])

    # 1.2. Add random noise to original matrix
    response <- AddNoise(response, method = "random")
    response.mat <- response %>%
      dplyr::select(conc_r, conc_c, response) %>%
      reshape2::acast(conc_r ~ conc_c, value.var = "response")

    # 2. Do calculation on matrix
    tmp <- CalculateMat(response.mat = response.mat)
    tmp <- lapply(tmp, function(x){
      x$block_id = rep(blocks[i], nrow(x))
      return(x)
      })

    # 3. Add information to summary table
    info <- response %>%
      dplyr::select(drug_row, drug_col, cell_line_name,
                    conc_r_unit, conc_c_unit) %>%
      unique()

    if (nrow(info) > 1) {
      warning("The summary data of block ", blocks[i], " contains ", nrow(info),
              " different versions.")
    } else if (nrow(info) < 1) {
      warning("The summary data of block ", blocks[i], " is missing.")
    }

    tmp$summary <- cbind.data.frame(info, tmp$summary)

    # 4. fill drug names to curve table

    tmp$curve$drug_col <- rep(NA, 2)
    tmp$curve$drug_col[which(tmp$curve$dim ==
                               "col")] <- as.character(info$drug_col)

    tmp$curve$drug_row <- rep(NA, 2)
    tmp$curve$drug_row[which(tmp$curve$dim ==
                               "row")] <- as.character(info$drug_row)

    tmp$curve <- dplyr::select(tmp$curve, block_id, drug_row, drug_col, b, c, d, e, model)

    # # 4. Append new tables to container
    # synergy <- rbind.data.frame(synergy, tmp$synergy)
    # surface <- rbind.data.frame(surface, tmp$surface)
    # curve <- rbind.data.frame(curve, tmp$curve)
    # summary <- rbind.data.frame(summary, tmp$summary)

    # 4. collect tables
    result$synergy <- tmp$synergy
    result$surface <- tmp$surface
    result$summary <- tmp$summary
    result$curve <- tmp$curve

    # Clean temporary file
    tmp <- list()
    info <- data.frame()
    response <- data.frame()
    respons.mat <- matrix()

    gc()

    result
  }

  res2 <- list()
  res2$synergy <- Reduce(function(x, y) {rbind.data.frame(x, y)}, lapply(res, "[[" , "synergy"))
  res2$surface <- Reduce(function(x, y) {rbind.data.frame(x, y)}, lapply(res, "[[" , "surface"))
  res2$summary <- Reduce(function(x, y) {rbind.data.frame(x, y)}, lapply(res, "[[" , "summary"))
  res2$curve <- Reduce(function(x, y) {rbind.data.frame(x, y)}, lapply(res, "[[" , "curve"))

  return(res2)
}
