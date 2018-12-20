# CSS - GenerateSummary
# - version 1: Joseph
# The following script defines the wrapper function generateSummary() which
# accepts a csv file as input (the whole dataset) and produces the a csv
# summary output file via a series of structure remodeling steps and
# sequential function calls
# ic50_row: the selected (IC50) concentration of the row drug
# ic50_col: the selected (IC50) concentration of the col drug
# dss_row: the dss for the row drug, i.e the first column
# dss_col: the dss for the col drug, i.e the first row
# css_row: a particular column selected according to ic50_col
# css_col: a particular row selected according to ic50_row

GenerateSummary <- function(bigdata) {
  # Select the needed columns for Steps 1&2
  neededcolumns1 <- c("block_id", "cell_line_name", "drug_row", "conc_r",
                      "drug_col", "conc_c", "response")
  raw = bigdata[, which(colnames(bigdata) %in% neededcolumns1)]

  #############################################################################
  # Step1: Create a summary IC50 table
  # calls createIC50Summary() which calls calculateIC50()
  print(noquote("Step 1"))
  print(noquote("Calculating IC50 concentrations..."))
  drugsIC50 <- createIC50Summary(raw)
  save(raw, drugsIC50, file = "step1.RData")

  print(noquote(""))
  print(noquote("An IC50 concentration has been assigned to each drug."))
  print(noquote(""))

  ##############################################################################
  # Step2: Summarize experiments by calculating four sensitivity scores for each
  # calls summarizeBlocks() which calls imputeData() and/or computeSensitivity()
  # which calls scoreCurve()
  print(noquote("Step 2"))
  print(noquote("Calculating sensitivity scores..."))
  sensscores = summarizeBlocks(raw, drugsIC50)
  save(raw, sensscores, file = "step2.RData")

  print(noquote(""))
  print(noquote("DSS and CSS values have been calculated per experiment."))
  print(noquote(""))

  ##############################################################################
  # Select the needed columns for Step 3 and drop the unwanted rows (those
  # containing single-agent data)
  neededcolumns2 <- c("block_id", "conc_r", "conc_c", "synergy_zip",
                     "synergy_hsa", "synergy_bliss", "synergy_loewe")
  syn <- bigdata[, which(colnames(bigdata) %in% neededcolumns2)]
  syn <- syn[which(syn$'conc_r' != 0 & syn$'conc_c' != 0),
                   c("block_id", "synergy_zip", "synergy_hsa",
                     "synergy_bliss", "synergy_loewe")]

  ##############################################################################
  # Step3: Summarize synergies per experiment
  # Calls summarizeSynergies()
  print(noquote("Step 3"))
  print(noquote("Summarizing synergy values..."))
  synscores <- summarizeSynergies(syn)
  save(bigdata, synscores, file = "step3.RData")
  print(noquote(""))
  print(noquote("Synergy scores have been summarized per experiment."))
  print(noquote(""))

  ##############################################################################
  # Step4: combine sensscores and synscores and design the output as required
  print(noquote("Step 4"))
  summarytable <- merge(sensscores,synscores, by = 0, all = TRUE)
  summarytable <- summarytable[order(as.numeric(summarytable$Row.names)), ]
  rownames(summarytable) <- NULL
  colnames(summarytable) <- c('block_id', 'cell_line_name', 'drug_row_name',
                              'drug_col_name', 'ic50_row', 'ic50_col',
                              'dss_row', 'dss_col', 'css_row', 'css_col', 'css',
                              'synergy_zip', 'synergy_hsa', 'synergy_bliss',
                              'synergy_loewe')
  # summarytable <- summarytable[,c(1,5,6,11,7:10,12,14,15,13,3,4,2)]
  summarytable <- summarytable[, c("block_id", "ic50_row", "ic50_col", "css",
                                   "dss_row", "dss_col", "css_row", "css_col",
                                   "synergy_zip", "synergy_bliss",
                                   "synergy_loewe", "synergy_hsa",
                                   "drug_row_name", "drug_col_name",
                                   "cell_line_name")]
  save(summarytable, file = "step4.RData")
  print(noquote("Final table has been successfully built."))

  ##############################################################################
  return(summarytable)
}
