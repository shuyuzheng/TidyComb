# ReshapeData
# Version 3: for analying all the scores simultaneously
ReshapeData3 = function (data) {
  if (!all(c("block_id", "drug_row", "drug_col", "row", "col", "response",
             "synergy_zip", "synergy_loewe", "synergy_bliss", "synergy_hsa",
             "conc_r", "conc_c", "conc_r_unit", "conc_c_unit") %in%
           colnames(data)))
    stop("The input data must contain the following columns:",
         " block_id, drug_row, drug_col, row, col, response, synergy_zip, ",
         "synergy_loewe, synergy_bliss, synergy_hsa \n",
         "conc_r, conc_c, conc_r_unit, conc_c_unit")
  id.drug.comb <- unique(data$block_id)
  dose.response.mats.response <- list()
  dose.response.mats.zip <- list()
  dose.response.mats.loewe <- list()
  dose.response.mats.bliss <- list()
  dose.response.mats.hsa <- list()

  drug.pairs <- data.frame(drug.row = character(length(id.drug.comb)),
                           drug.col = character(length(id.drug.comb)),
                           concRUnit = character(length(id.drug.comb)),
                           concCUnit = character(length(id.drug.comb)),
                           blockIDs = numeric(length(id.drug.comb)),
                           stringsAsFactors = FALSE)
  for (i in 1:length(id.drug.comb)) {
    # cat('\r', i)
    tmp.mat <- data[which(data$block_id == id.drug.comb[i]),]

    conc.col <- tmp.mat$conc_c[which(tmp.mat$row == 1)]
    conc.col <- conc.col[order(tmp.mat$col[which(tmp.mat$row == 1)])]
    conc.row <- tmp.mat$conc_r[which(tmp.mat$col == 1)]
    conc.row <- conc.row[order(tmp.mat$row[which(tmp.mat$col == 1)])]
    # replicates are averaged
    response.mat.response <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                             value.var = "response",
                                             function(x) mean(x, na.rm = T))

    response.mat.zip <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                        value.var = "synergy_zip",
                                        function(x) mean(x, na.rm = TRUE))
    response.mat.loewe <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                          value.var = "synergy_loewe",
                                          function(x) mean(x, na.rm = TRUE))
    response.mat.bliss <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                          value.var = "synergy_bliss",
                                          function(x) mean(x, na.rm = TRUE))
    response.mat.hsa <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                        value.var = "synergy_hsa",
                                        function(x) mean(x, na.rm = TRUE))

    conc.runit <- unique(tmp.mat$conc_r_unit)
    conc.cunit <- unique(tmp.mat$conc_c_unit)

    drug.row <- unique(tmp.mat$drug_row)
    drug.col <- unique(tmp.mat$drug_col)

    drug.pairs$drug.row[i] <- drug.row
    drug.pairs$drug.col[i] <- drug.col
    drug.pairs$concRUnit[i] <- conc.runit
    drug.pairs$concCUnit[i] <- conc.cunit
    dose.response.mats.response[[i]] <- response.mat.response
    dose.response.mats.zip[[i]] <- response.mat.zip
    dose.response.mats.loewe[[i]] <- response.mat.loewe
    dose.response.mats.bliss[[i]] <- response.mat.bliss
    dose.response.mats.hsa[[i]] <- response.mat.hsa
  }
  drug.pairs$blockIDs <- id.drug.comb
  return(list(dose.response.mats.response = dose.response.mats.response,
              dose.response.mats.zip = dose.response.mats.zip,
              dose.response.mats.loewe = dose.response.mats.loewe,
              dose.response.mats.bliss = dose.response.mats.bliss,
              dose.response.mats.hsa = dose.response.mats.hsa,
              drug.pairs = drug.pairs))
}
