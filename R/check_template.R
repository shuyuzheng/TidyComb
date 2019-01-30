CheckTemplate <- function(template) {

  # Checking necessary columns for calculation
  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r",
             "conc_c", "conc_r_unit", "conc_c_unit","cell_line_name") %in%
           colnames(response))) {
    stop("The template must contain the following columns: ",
         "block_id, response, conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
         "cell_line_name, drug_row, drug_col.")
  }

  # Checking NA rows or NA columns
  missing.response <- template %>%
    select(block_id, drug_col, )


}