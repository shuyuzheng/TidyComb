# CSS - GenerateSummary
# Version 2 - Jing
# The following script defines the wrapper function generateSummary() which
# accepts a csv file as input (the whole dataset) and produces the a csv
# summary output file via a series of structure remodeling steps and sequential
# function calls
# ic50_row: the selected (IC50) concentration of the row drug
# ic50_col: the selected (IC50) concentration of the col drug
# dss_row: the dss for the row drug, i.e the first column
# dss_col: the dss for the col drug, i.e the first row
# css_row: a particular column selected according to ic50_col
# css_col: a particular row selected according to ic50_row

GenerateSummary2 <- function(response) {
  sensscores <- unique(response[, c("block_id", "cell_line_id",
                                    "drug_row_id", "drug_col_id")])
  rownames(sensscores) <- sensscores$block_id
  m <- sensscores$block_id
  len<- length(m)

  sensscores <- sensscores[, -1]
  sensscores$IC_row  <- NA
  sensscores$IC_col  <- NA
  sensscores$DSS_row <- NA
  sensscores$DSS_col <- NA
  sensscores$CSS_row <- NA
  sensscores$CSS_col <- NA

  options(show.error.messages = FALSE)
  # each block ID defines an experiment, i.e., a combination for a particular
  # cell line
  for (i in 1:len) {
    # cat("\r", i)
    print(m[i])
    # extract comb info for this block
    tempmaindf <- response[which(response$'block_id' == m[i]), ]
    # assume no replicates
    response.mat <- reshape2::acast(tempmaindf, conc_r ~ conc_c,
                                    value.var = "response")

    # missing value imputation - NB! the imputation only stay within the
    # function find the average of the neighboring cells in the matrix
    x <- array(c(rbind(response.mat[-1,], NA),
                rbind(NA, response.mat[-nrow(response.mat),]),
                cbind(response.mat[,-1], NA),
                cbind(NA, response.mat[, -ncol(response.mat)])),
              dim=c(nrow(response.mat),ncol(response.mat),4))
    x.imp <- apply(x, c(1, 2),function(x) mean(x, na.rm = TRUE))
    index.na <- is.na(response.mat)
    response.mat[index.na] <- x.imp[index.na]

    # one more round
    x <- array(c(rbind(response.mat[-1,], NA),
                rbind(NA, response.mat[-nrow(response.mat),]),
                cbind(response.mat[,-1], NA),
                cbind(NA, response.mat[, -ncol(response.mat)])),
              dim=c(nrow(response.mat),ncol(response.mat),4))
    x.imp <- apply(x, c(1, 2), function(x) mean(x, na.rm = TRUE))
    index.na <- is.na(response.mat)
    response.mat[index.na] <- x.imp[index.na]

    scale <- 10^-10
    noise <- row(response.mat) * col(response.mat) * scale
    response.mat <- response.mat + noise

    # DSS_c
    tempdf_c <- data.frame(conc = as.numeric(colnames(response.mat)),
                           response = response.mat[1, ]) # First row
    DSS_c <- computeSensitivity(tempdf_c)

    # DSS_r
    tempdf_r <- data.frame(conc = as.numeric(rownames(response.mat)),
                           response = response.mat[, 1]) # First col
    DSS_r <- computeSensitivity(tempdf_r)

    # Relative IC50
    ic_c <- calculateIC50(tempdf_c) # ic50 for col drug
    ic_r <- calculateIC50(tempdf_r) # ic50 for row drug

    # CSS
    res <- imputeData2(response.mat, ic_c, ic_r)
    # a particular row selected according to ic50_row
    CSS_c <- computeSensitivity(res$tempcf_c)
    # a particular column selected according to ic50_col
    CSS_r <- computeSensitivity(res$tempcf_r)

    # assign six measurements to sensscores
    sensscores[i, c("IC_row", "IC_col", "DSS_row",
                    "DSS_col", "CSS_row", "CSS_col")] <- c(ic_r, ic_c, DSS_r,
                                                           DSS_c, CSS_r, CSS_c)

  }
  options(show.error.messages = TRUE)

  # save(response, sensscores, file = "step1.RData")
  sensscores$CSS <- NA
  for(i in 1:nrow(sensscores)) {
    sensscores$CSS[i] <- mean(c(sensscores$CSS_r[i], sensscores$CSS_c[i]),
                              na.rm = TRUE)
  }

  neededcolumns2 <- c("block_id", "conc_r", "conc_c", "synergy_zip",
                      "synergy_hsa", "synergy_bliss", "synergy_loewe")
  syn <- response[, which(colnames(response) %in% neededcolumns2)]
  syn <- syn[which(syn$'conc_r'!= 0 & syn$'conc_c'!= 0),
          c("block_id", "synergy_zip", "synergy_hsa",
            "synergy_bliss", "synergy_loewe")]

  # print(noquote("Summarizing synergy values..."))
  synscores <- summarizeSynergies(syn)
  # save(response, synscores, file = "step2.RData")
  # print(noquote(""))
  # print(noquote("Synergy scores have been summarized per experiment."))
  # print(noquote(""))

  ##############################################################################
  # Step3: combine sensscores and synscores and design the output as required
  # print(noquote("Step 3"))
  summarytable <- merge(sensscores, synscores, by = 0, all = TRUE)
  summarytable <- summarytable[order(as.numeric(summarytable$Row.names)), ]
  rownames(summarytable) <- NULL
  colnames(summarytable) <- c('block_id', 'cell_line_id', 'drug_row_id',
                              'drug_col_id', 'ic50_row', 'ic50_col',
                              'dss_row', 'dss_col', 'css_row', 'css_col', 'css',
                              'synergy_zip', 'synergy_hsa', 'synergy_bliss',
                              'synergy_loewe')
  # summarytable <- summarytable[,c(1,5,6,11,7:10,12,14,15,13,3,4,2)]
  summarytable <- summarytable[, c("block_id", "ic50_row", "ic50_col", "css",
                                   "dss_row", "dss_col", "css_row", "css_col",
                                   "synergy_zip", "synergy_bliss",
                                   "synergy_loewe", "synergy_hsa",
                                   "drug_row_id", "drug_col_id",
                                   "cell_line_id")]
  summarytable$conc_r_unit <- unique(response$conc_r_unit)
  summarytable$conc_c_unit <- unique(response$conc_c_unit)
  # save(summarytable, file = "step3.RData")
  # print(noquote("Final table has been successfully built."))

  ##############################################################################
  return(summarytable)
}
