# Function to generate the surface table

# response_with_scores must contain the following columns
# c("block_id", "conc_r","conc_c", "response", "drug_row", "drug_col",
#   "conc_r_unit", "conc_c_unit", "cell_line_name", "row", "col", "synergy_zip",
#   "synergy_hsa", "synergy_bliss", "synergy_loewe")

# Version 3: Analyze all the scores simultanesouly
GenerateSurface3 <- function(response_with_scores) {
  if (!all(c("block_id", "drug_row", "drug_col", "row", "col", "response",
             "synergy_zip", "synergy_loewe", "synergy_bliss", "synergy_hsa",
             "conc_r", "conc_c", "conc_r_unit", "conc_c_unit") %in%
           colnames(response_with_scores)))
    stop("The input data must contain the following columns:",
         " block_id, drug_row, drug_col, row, col, response, synergy_zip, ",
         "synergy_loewe, synergy_bliss, synergy_hsa \n",
         "conc_r, conc_c, conc_r_unit, conc_c_unit")

  m <- unique(response_with_scores$block_id)
  len <- length(m)

  data2 <- ReshapeData3(response_with_scores)
  surface <- list()

  options(show.error.messages = FALSE)
  for(i in 1:len) {
    # cat('\r', i)
    print(m[i])
    # flush.console()
    index <- which(data2$drug.pairs$blockIDs == m[i])

    # surface krigging
    mat.response <- smoothing(data2$dose.response.mats.response[[index]])
    mat.zip <- smoothing(data2$dose.response.mats.zip[[index]])
    mat.loewe <- smoothing(data2$dose.response.mats.loewe[[index]])
    mat.bliss <- smoothing(data2$dose.response.mats.bliss[[index]])
    mat.hsa <- smoothing(data2$dose.response.mats.hsa[[index]])

    # list format
    mat2.response <- stats::setNames(reshape2::melt(mat.response),
                                     c('conc_r', 'conc_c', 'response'))
    mat2.zip <- stats::setNames(reshape2::melt(mat.zip),
                               c('conc_r', 'conc_c', 'synergy_zip'))
    mat2.loewe <- stats::setNames(reshape2::melt(mat.loewe),
                                 c('conc_r', 'conc_c', 'synergy_loewe'))
    mat2.bliss <- stats::setNames(reshape2::melt(mat.bliss),
                                 c('conc_r', 'conc_c', 'synergy_bliss'))
    mat2.hsa <- stats::setNames(reshape2::melt(mat.hsa),
                                c('conc_r', 'conc_c', 'synergy_hsa'))

    mat2 <- cbind(mat2.response,
                  synergy_zip = mat2.zip$synergy_zip,
                  synergy_loewe = mat2.loewe$synergy_loewe,
                  synergy_bliss = mat2.bliss$synergy_bliss,
                  synergy_hsa = mat2.hsa$synergy_hsa)

    mat2$block_id <- data2$drug.pairs$blockIDs[index]
    mat2$drug_row <- data2$drug.pairs$drug.row[index]
    mat2$drug_col <- data2$drug.pairs$drug.col[index]
    mat2$conc_r_unit <- data2$drug.pairs$concRUnit[index]
    mat2$conc_c_unit <- data2$drug.pairs$concCUnit[index]


    surface[[i]] <- mat2
  }
  options(show.error.messages = TRUE)
  surface_table <- do.call(rbind, surface)
  return(surface_table)
}
