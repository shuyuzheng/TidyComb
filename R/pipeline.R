# Function to take the uploaded response data
# pipeline(response)


#' The pipeline for calculation
#'
#' @param response A data frame. It must contain columns: "block_id",
#' "drug_row_id", "drug_col_id", "response", "conc_r", "conc_c", "conc_r_unit",
#'  "conc_c_unit","cell_line_id"
#'
#' @return A list contains 6 data frames: "response_with_scores", "response",
#' "summary", "surface", "curve".
#' @export
#'
pipeline <- function(response) {

  options(scipen = 999) # suppress the scientific notation

  if (!all(c("block_id", "drug_row_id", "drug_col_id", "response", "conc_r",
             "conc_c", "conc_r_unit", "conc_c_unit","cell_line_id") %in%
           colnames(response)))
    stop("The input data must contain the following columns: ",
         "block_id, response, conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
         "cell_line_id, drug_row_id, drug_col_id.")

  print(noquote('Generating scores...'))
  response_with_scores <- GenerateScore(response)
  response_out <- response_with_scores[, c("block_id", "conc_r", "conc_c",
                                           "response", "synergy_zip",
                                           "synergy_bliss", "synergy_loewe",
                                           "synergy_hsa")]

  print(noquote('Generating summary...'))
  summary <- GenerateSummary2(response_with_scores)

  print(noquote('Generating surface...'))
  surface <- GenerateSurface3(response_with_scores)

  print(noquote('Generating curve...'))
  curve <- GenerateCurve(response_with_scores)

  print(noquote('Ready.'))
  return(list(response_with_scores = response_with_scores,
              response = response_out,
              summary = summary,
              surface = surface,
              curve = curve))
}