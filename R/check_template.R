# TidyComb
# Function for checking input template
# Copyright: Shuyu Zheng

#' Check the variables of "template" file
#'
#' Function \code{CheckTemplate} checks the variables in input \code{template}
#' to see whether all required elements have been provided.
#'
#' The \code{template} file must have following variables:
#' \itemize{
#'   \item \strong{block_id} Identifier for each drug combination reaction
#'   block.
#'   \item \strong{drug_row} Identifier for drugs added to rows in certain
#'   block.
#'   \item \strong{drug_col} Identifier for drugs added to columns in certain
#'   block.
#'   \item \strong{response} Response value (inhibition rate for cell growth or
#'   other
#'   type of signals) detected from certain well in test blocks.
#'   \item \strong{conc_r} Concentration of drugs added to rows in certain
#'   block.
#'   \item \strong{conc_c} Concentration of drugs added to columns in certain
#'   block.
#'   \item \strong{conc_r_unit} The unit used to measure concentration of drugs
#'   added in rows.
#'   \item \strong{conc_c_unit} The unit used to measure concentration of drugs
#'   added in columns.
#'   \item \strong{cell_line_name} The name of cell line in which drug
#'   combination screen was performed.
#' }
#'
#' @param template The template table which is provided by user for uploading
#' data to DrugComb database.
#'
#' @return If the template lack some of required elements, the program will be
#' shut down with an error. If the the template is all right, the program will
#' continue with a message tell user the template is OK.
#'
#' @export
CheckTemplate <- function(template) {

  # Checking necessary columns for calculation
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r",
             "conc_c", "conc_r_unit", "conc_c_unit","cell_line_name") %in%
           colnames(response))) {
    stop("The template must contain the following columns: ",
         "block_id, response, conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
         "cell_line_name, drug_row, drug_col.")
  } else {
    message("The template is OK and the data process continues.")
  }
}