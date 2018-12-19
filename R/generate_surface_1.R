# Function to generate the surface table

# response_with_scores must contain the following columns
# c("block_id","conc_r","conc_c", "response", "drug_row", "drug_col",
#   "conc_r_unit", "conc_c_unit", "cell_line_name", "row", "col", "synergy_zip",
#   "synergy_hsa", "synergy_bliss", "synergy_loewe")

GenerateSurface = function(response_with_scores, model){
  #model %in% c("response", "synergy_loewe", "synergy_bliss", "synergy_hsa",
  #             "synergy_zip")

  if (!model %in% c("response", "synergy_loewe", "synergy_bliss", "synergy_hsa",
                    "synergy_zip")) {
    stop("The model can only be one of the following:",
         " response, synergy_loewe, synergy_bliss, synergy_hsa, synergy_zip")
  }
  m <- unique(response_with_scores$block_id)
  len <- length(m)

  data <- cbind(response_with_scores[,
                                  which(names(response_with_scores) == model)],
                response_with_scores[, c("block_id", "drug_row", "drug_col",
                                         "row", "col", "conc_r", "conc_c",
                                         "conc_r_unit", "conc_c_unit")])
  names(data)[1] <- "response"

  data2 <- ReshapeData2(data, data.type = "inhibition")
  surface <- list()
  # options(show.error.messages = F)
  for(i in 1:len){
    cat('\r', i)
    # flush.console()
    index <- which(data2$drug.pairs$blockIDs == m[i])

    # surface krigging
    mat1 <- tryCatch ({ # in case error happens, return the score of NA
        mat1 <- smoothing(data2$dose.response.mats[[index]])
        # my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
        # heatmap.2(mat1, Rowv = F, Colv = F, dendrogram = 'none',
        #           trace = 'none', col = my_palette, density.info = "none")
      }, error = function(cond) {
        print(i);
        mat1 <- data2$dose.response.mats[[index]]
        mat1[] <- 0
        mat1 <- smoothing(mat1)
        return(mat1)
      }
    )

    # list format
    mat2 <- stats::setNames(melt(mat1), c('conc_r', 'conc_c', 'response'))
    mat2$block_id <- data2$drug.pairs$blockIDs[index]
    mat2$drug_row <- data2$drug.pairs$drug.row[index]
    mat2$drug_col <- data2$drug.pairs$drug.col[index]
    mat2$conc_r_unit <- data2$drug.pairs$concRUnit[index]
    mat2$conc_c_unit <- data2$drug.pairs$concCUnit[index]
    surface[[i]]<- mat2
  }
  # options(show.error.messages = TRUE)
  surface_table <- do.call(rbind, surface)
  names(surface_table)[which(names(surface_table)=="response")] <- model
  return(surface_table)
}
