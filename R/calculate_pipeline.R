# TidyComb
# Functions for pipe all calculation togather
#
# Functions in this page:
#
# CalculateMat: Calculate scores from one dose-response matrix.
# CalculateTemplate: Calculate scores from template file which might contain
#                    several drug combination blocks.
# CalculateTemplateDebug: Debug version of CalculateTemplate.
# ParCalculateTemplate: Parallel calculate scores from template file. It using
#                    multi-core method to do parallel. (It might not work under
#                    Windows system)

#' Chain all calculation about one dose-response matrix
#'
#' \code{CalculateMat} chains all calculations about one dose-response matrix (
#' one drug-drug interaction block) together. The calculations includes:
#' dose-response curve fitting, synergy scores (ZIP, Bliss, Loewe,
#' HSA, S), drug sensitivity (RI, CSS), generate drug-drug response surface and
#' generating summary scores for each block.
#'
#' The steps for calculation:
#' \enumerate{
#'   \item Pre-process Matrix
#'     \enumerate{
#'       \item Impute for missing values (with the average of values from the
#'           nearest four cells) in original matrix by using function
#'           \code{\link[synergyfinder]{ImputeNA}}.
#'       \item Add noise(A small random number ranging from 0 to 0.001) to
#'           original matrix by using function \code{\link[synergyfinder]{AddNoise}}.
#'       )
#'       \item Correct baseline using function
#'       \code{\link[synergyfinder]{CorrectBaseLine}} with the method
#'       selected by parameter \code{correction}.
#'     }
#'   \item Single drug process
#'     \enumerate{
#'       \item Extract and fitting single drugs.
#'       \item Extract coeficients from fitted model. (b, c, d, e, IC50)
#'       \item Calculate RI(Relative inhibition for single drug) with function
#'       \code{CalculateSens}
#'     }
#'   \item Whole response matrix process
#'     \enumerate{
#'       \item Calculate Synergy Scores with function \code{\link[synergyfinder]{ZIP}},
#'       \code{\link[synergyfinder]{Bliss}}, \code{\link[synergyfinder]{HSA}},
#'       \code{\link[synergyfinder]{Loewe}} in \code{synergyfinder} package.
#'       \item Calculate Surface(The landscape of response, synergy scores)
#'       \item Calculate CSS(drug combination sensitivity score), S(synergy
#'       score calculated from CSS and IR)
#'     }
#'   \item Summarize and generate surface
#' }
#' @param response.mat A matrix contain the drug combination reaponse value.
#' Its column names are doses of drug added along columns. Its row name are
#' doses of drug added along rows. \cr
#' \strong{Note}: the matrix should be sorted by: 1. The concentrations along
#' the column increase \emph{from left to right}; 2. The concentrations along
#' the row increase \emph{from top to bottom}.
#'
#' @param noise a logical value. It indicates whether or not adding noise to
#' to the "inhibition" values in the matrix. Default is \code{TRUE}.
#'
#' @param correction a string. It indicates which method used by function
#' \code{\link[synergyfinder]{CorrectBaseLine}} for base line correction.
#'   \itemize{
#'     \item \code{non}: no baseline correction;
#'     \item \code{par}: only correct base line on negative values in the matrix;
#'     \item \code{all}: correct base line on all the values in the matrix.
#'   }
#'
#' @param summary.only a logical value. If it is \code{TRUE} then only summary
#' table is calculated and returned, otherwise, for tables will be return.
#' Default setting is \code{FALSE}.
#'
#' @return A list contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified inhibition value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, ri, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed inhibition values and
#'     synergy scores for plotting the scores' landscape.
#'  }
#'  If \code{summary.only} is \code{TRUE}, it will return only the "summary"
#'  data frame.
#'
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data <- read.csv(system.file("template.csv", package = "TidyComb"),
#'                  stringsAsFactors = FALSE)
#' response.mat <- reshape2::acast(conc_r~conc_c, value.var = "inhibition",
#'                                 data = data[data$block_id == 1, ])
#' res <- CalculateMat(response.mat)
CalculateMat <- function(response.mat, noise = TRUE, correction = "non",
                         summary.only = FALSE) {

  options(scipen = 999)

  # 0. Check input data
  if (!is.matrix(response.mat)) {
    warning("Input data is not 'matrix'. Trying to transform it.")
    response.mat <- as.matrix(response.mat)
  }

  # 1. Pre-processing
  # 1.1 Impute data until there is no missing value in the response matrix
  while (sum(is.na(response.mat))) {
    response.mat <- synergyfinder::ImputeNA(response.mat)
  }

  # 1.2. Add random noise to original matrix
  set.seed(1)
  if (noise) {
    response.mat <- synergyfinder::AddNoise(response.mat)
  }

  # 1.3. Correct baseline with corresponding "method". Available methods are
  #      "non", "part", "all".
    response.mat <- synergyfinder::CorrectBaseLine(response.mat,
                                                   method = correction)
  # 2. Single drug process
  # 2.1. Fit single drug dose-response curve
  # drug_col
  drug.col <- synergyfinder::ExtractSingleDrug(response.mat, dim = "col")
  col.model <- synergyfinder::FitDoseResponse(drug.col)
  col.type <- as.character(synergyfinder::FindModelType(col.model))

  # drug_row
  drug.row <- synergyfinder::ExtractSingleDrug(response.mat, dim = "row")
  row.model <- synergyfinder::FitDoseResponse(drug.row)
  row.type <- as.character(synergyfinder::FindModelType(row.model))

  # 2.2 Extract coeficients
  # drug_col
  col.coe <- stats::coef(col.model)
  # drug_rowr
  row.coe <- stats::coef(row.model)

  if (!summary.only) {
    curves <- as.data.frame(rbind(col.coe, row.coe))
    colnames(curves) <- sub(":(Intercept)", "", colnames(curves), fixed = TRUE)
    curves$model <- c(col.type, row.type)
    curves$dim <- c("col", "row")
  }


  # 2.3 Calculate RI (using single drug response but without that at 0
  # concentration)
  col.ri <- CalculateSens(drug.col)
  row.ri <- CalculateSens(drug.row)

  # 3. whole matrix process
  # 3.1 Calculate synergyscores
  zip <- synergyfinder::ZIP(response.mat, drug.row.model = row.model,
                                     drug.col.model = col.model)
  loewe <- synergyfinder::Loewe(response.mat, drug.row.model = row.model,
                                drug.col.model = col.model)
  hsa <- synergyfinder::HSA(response.mat)
  bliss <- synergyfinder::Bliss(response.mat)

  synergy <- lapply(list(response.mat, zip, loewe, hsa, bliss), reshape2::melt)
  synergy <- Reduce(function(x, y) {
    merge(x = x, y = y, by = c("Var1", "Var2"))}, synergy)

  colnames(synergy) <- c("conc_r", "conc_c", "inhibition", "synergy_zip",
                         "synergy_loewe", "synergy_hsa", "synergy_bliss")

  if (!summary.only) {
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

    colnames(surface) <- c("conc_r", "conc_c", "inhibition", "synergy_zip",
                           "synergy_loewe", "synergy_hsa", "synergy_bliss")
  }


  # 3.3 Calculate CSS
  col.ic50 <- CalculateIC50(col.coe, col.type, max(drug.col$dose))
  row.ic50 <- CalculateIC50(row.coe, row.type, max(drug.row$dose))

  imputed.ic50 <- ImputeIC50(response.mat, col.ic50 = col.ic50,
                             row.ic50 = row.ic50)
  # a particular row selected according to ic50_row
  col.css <- CalculateSens(imputed.ic50$tempcf_c)
  # a particular column selected according to ic50_col
  row.css <- CalculateSens(imputed.ic50$tempcf_r)

  css <- mean(c(col.css, row.css), na.rm = TRUE)

  S <- css - sum(col.ri, row.ri)

  sum <- synergy %>%
    dplyr::filter(conc_r != 0 & conc_c != 0) %>%
    dplyr::select(-conc_r, -conc_c, -inhibition) %>%
    apply(2, mean, na.rm = TRUE)

  sum <- data.frame(t(sum), ic50_row = row.ic50 , ic50_col = col.ic50,
                    ri_row = row.ri, ri_col = col.ri, css_row = row.css,
                    css_col = col.css, css = css, S = S)
  if (summary.only) {
    return(sum)
  } else{
    res <- list(synergy = synergy, surface = surface,
                summary = sum, curve = curves)
    return(res)
  }

  # clean up
  gc()
  rm(.Random.seed)
}

#' Calculate Drug Combination data in template format
#'
#' @param template a dataframe which must contains following columns:
#'  \itemize{
#'   \item \emph{block_id}: (integer) the identifier for a drug combination.
#'    If mul-tiple drug combinations are present, e.g. in the standard 384-well
#'    platewhere 6 drug combinations are fitted, then the identifiers for each
#'    of themmust be unique.
#'    \item \emph{drug_col}: (character) the name of the drug on the columns in
#'    adose-response matrix.
#'    \item \emph{drug_row}: (character) the name of the drug on the rows in
#'    adose-response matrix.
#'    \item \emph{conc_c}: (numeric) the concentrations of the column drugs
#'    in a combination.
#'    \item \emph{conc_r}: (numeric) the concentrations of the row drugs
#'    in a combination.
#'    \item \emph{conc_c_unit}: (character) the unit of concentrations of the
#'    column drugs. It is typically nM or \eqn{\mu}M.
#'    \item \emph{conc_r_unit}: (character) the unit of concentrations of the
#'    row drugs. It is typically nM or \eqn{\mu}M.
#'    \item \emph{response}: (numeric) the effect of drug combinations at the
#'    concentra-tions specified by conc_r and conc_c. The effect must be
#'    normalized to \%inhibition based on the positive and
#'    negative controls. For a well-controlled experiment, the range of the
#'    response values is expected from 0 to 100. However, missing values or
#'    extreme values are allowed.
#'    \item \emph{cell_line_name}: (character) the name of cell lines on which
#'    the drug combination was tested.
#'  }
#'
#' @param summary.only a logical value. If it is \code{TRUE} then only summary
#' table is calculated and returned, otherwise, for tables will be return.
#' Default setting is \code{FALSE}.
#'
#' @param ... Other arguments required by nested functions. Some important
#' arguments are:
#'  \itemize{
#'    \item \code{correction} and \code{noise} inherited from function
#'          \code{\link{CalculateMat}};
#'    \item \code{Emin} and \code{Emax} inherited from function
#'          \code{\link{FitDoseResponse}}.
#' }
#'
#' @return A list. It contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified response value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, ri, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed response value and
#'     synergy scores of each drug dose response pair.
#'  }
#'
#' @author Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' data <- read.csv(system.file("template.csv", package = "TidyComb"),
#'                  stringsAsFactors = FALSE)
#' res <- CalculateTemplate(data)
CalculateTemplate <- function(template, summary.only=FALSE) {
  CheckTemplate(template)

  blocks <- unique(template$block_id)

  # generate container
  summary <- data.frame(block_id = integer(), synergy_zip = numeric(),
                        synergy_bliss = numeric(), synergy_hsa = numeric(),
                        synergy_loewe = numeric(), ic50_row = numeric() ,
                        ic50_col = numeric(), ri_row = numeric(),
                        ri_col = numeric(), css_row = numeric(),
                        css_col = numeric(), css = numeric(), S = numeric(),
                        stringsAsFactors = FALSE)

  if (!summary.only) {
    synergy <- data.frame(block_id = integer(), conc_r = numeric(),
                          conc_c = numeric(), inhibition = numeric(),
                          synergy_zip = numeric(), synergy_bliss = numeric(),
                          synergy_loewe = numeric(), synergy_hsa = numeric(),
                          stringsAsFactors = FALSE)
    surface <- synergy
    curve <- data.frame(block_id = integer(), b = numeric(),
                        c = numeric(), d = numeric(), e = numeric(),
                        model = numeric(), drug.row = numeric(),
                        drug.col = numeric(), stringsAsFactors = FALSE)
  }

  for (block in blocks) {
    # 1. Generate response matrix for each block
    response <- template %>%
      dplyr::filter(block_id == block)

    response.mat <- response %>%
      dplyr::select(conc_r, conc_c, inhibition) %>%
      reshape2::acast(conc_r ~ conc_c, value.var = "inhibition")

    # 2. Do calculation on matrix (with error control)
    tmp <- tryCatch({
      CalculateMat(response.mat = response.mat,
                   summary.only = summary.only)
    }, error = function(e) {
      print(block)
      traceback()
    })

    if (summary.only) {
      tmp$block_id <- rep(block, nrow(tmp))
    } else {
      tmp <- lapply(tmp, function(x){
        x$block_id = rep(block, nrow(x))
        return(x)
      })
    }


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

    if (!summary.only) {
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
    } else {
      tmp <- cbind.data.frame(info, tmp)
      summary <- rbind.data.frame(summary, tmp)
    }



    # Clean temporary file
    tmp <- NULL
    info <- data.frame()
    response <- data.frame()
    response.mat <- matrix()
  }

  if (summary.only) {
    return(summary)
  } else {
    curve <- dplyr::select(curve, block_id, drug_row, drug_col, b, c, d, e, model)

    return(list(synergy = synergy, surface = surface,
                curve = curve, summary = summary))
  }
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
#' "drug_row", "drug_col", "inhibition", "conc_r", "conc_c", "conc_r_unit",
#' "conc_c_unit","cell_line_name", "drug_row", "drug_col" are reqired.
#'
#' @param cores A integer. It indicates number of cores would be allocated to
#' the parallel processed
#'
#' @param summary.only a logical value. If it is \code{TRUE} then only summary
#' table is calculated and returned, otherwise, for tables will be return.
#' Default setting is \code{FALSE}.
#'
#' @param ... Other arguments required by nested functions. Some important
#' arguments are:
#'  \itemize{
#'    \item \code{impute} and \code{noise} inherited from function
#'          \code{CalculateMat};
#'    \item \code{method} inherited from function \code{CorrectBaseLine};
#'    \item \code{Emin} and \code{Emax} inherited from function
#'          \code{FitDoseResponse}.
#' }
#'
#' @return A list. It contains 4 tables:
#'   \itemize{
#'     \item \strong{synergy} It contains the modified inhibition value and 4
#'     type of synergy scores of each drug dose response pair.
#'     \item \strong{summary} It contains summarized information of each
#'     blocks: synergy scores, css, ri, S
#'     \item \strong{curve} It contains the coefficients from single drug dose
#'     response curve.
#'     \item \strong{surface} It contains the smoothed inhibition value and
#'     synergy scores of each drug dose response pair.
#'  }
#'
#' @author Shuyu Zheng{shuyu.zheng@helsinki.fi}
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#'
#' @export
#' @examples
#' data <- read.csv(system.file("template.csv", package = "TidyComb"),
#'                  stringsAsFactors = FALSE)
#' res <- ParCalculateTemplate(data)
ParCalculateTemplate <- function(template, cores = 1, summary.only = FALSE) {
  template <- CheckTemplate(template)
  blocks <- unique(template$block_id)

  cl <- parallel::makeForkCluster(cores)
  doParallel::registerDoParallel(cl)

  res <- foreach::foreach (i = 1:length(blocks)) %dopar% {
    set.seed(1)
    # 1. Generate response matrix for each block
    result <- multiResultClass()
    response <- template %>%
      dplyr::filter(block_id == blocks[i])

    response.mat <- response %>%
      dplyr::select(conc_r, conc_c, inhibition) %>%
      reshape2::acast(conc_r ~ conc_c, value.var = "inhibition")

    # 2. Do calculation on matrix
    tmp <- CalculateMat(response.mat = response.mat,
                        summary.only = summary.only)
    if (summary.only) {
      tmp$block_id = rep(blocks[i], nrow(tmp))
    } else {
      tmp <- lapply(tmp, function(x){
        x$block_id = rep(blocks[i], nrow(x))
        return(x)
      })
    }

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

    if (summary.only) {
      result <- cbind.data.frame(info, tmp)
    } else {
      tmp$summary <- cbind.data.frame(info, tmp$summary)

      # 4. fill drug names to curve table

      tmp$curve$drug_col <- rep(NA, 2)
      tmp$curve$drug_col[which(tmp$curve$dim ==
                                 "col")] <- as.character(info$drug_col)

      tmp$curve$drug_row <- rep(NA, 2)

      tmp$curve$drug_row[which(tmp$curve$dim ==
                                 "row")] <- as.character(info$drug_row)

      tmp$curve <- dplyr::select(tmp$curve, block_id, drug_row, drug_col,
                                 b, c, d, e, model)

      # 5. collect tables
      result$synergy <- tmp$synergy
      result$surface <- tmp$surface
      result$summary <- tmp$summary
      result$curve <- tmp$curve
    }

    # Clean temporary file
    tmp <- NULL
    info <- data.frame()
    response <- data.frame()
    response.mat <- matrix()
    rm(.Random.seed)
    return(result)
  }

  # Stop the clusters
  parallel::stopCluster(cl)

  # Collecting the results
  if (summary.only) {
    res2 <- Reduce(function(x, y) {rbind.data.frame(x, y)}, res)
  } else {
    res2 <- list()
    res2$synergy <- Reduce(function(x, y) {rbind.data.frame(x, y)},
                           lapply(res, "[[" , "synergy"))
    res2$surface <- Reduce(function(x, y) {rbind.data.frame(x, y)},
                           lapply(res, "[[" , "surface"))
    res2$summary <- Reduce(function(x, y) {rbind.data.frame(x, y)},
                           lapply(res, "[[" , "summary"))
    res2$curve <- Reduce(function(x, y) {rbind.data.frame(x, y)},
                         lapply(res, "[[" , "curve"))
  }

  return(res2)
}
