CheckTemplate <- function(template){
  # Column names
  missing_col <- setdiff(c("block_id", "drug_row", "drug_col", "inhibition",
                           "conc_r", "conc_c", "conc_r_unit", "conc_c_unit",
                           "cell_line_name", "drug_row", "drug_col"),
                         colnames(template))
  if (length(missing_col) != 0)
    stop('Column(s) "', paste(missing_col, collapse = '", "'),
         '" is(are) missing from uploaded data. Please check input file and',
         "re-upload the data.")

  # Duplicate combinations in one block
  block_dup <- duplicated(template[,c("conc_c", "conc_r", "block_id")])
  if (sum(block_dup) > 0) {
    stop("There are duplicated ('conc_r'-'conc_c') combinations in the block: ",
         paste(unique(template$block_id[block_dup]), collapse = ", "),
         ". Analysis function can not handle replicate data right now. Please",
         "remove the duplicated rows and re-upload the data.")
  }
}