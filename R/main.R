# Function to take the uploaded response data
# pipeline(response)


pipeline <- function(response) {
  # library(openxlsx)
  # library(plyr)
  # library(reshape2)
  # library(drc)
  # library(nleqslv)
  # library(gplots)
  # library(dplyr)
  # library(SpatialExtremes)
  #
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_log.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_log2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData3.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/FittingSingleDrug2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateSynergy2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ZIP2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Loewe2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Bliss2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/HSA2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/BaselineCorrectionSD2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/d1.fun.R') # used in loewe # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/d2.fun.R') # used in loewe # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/eq.R') # used in loewe # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateScore.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface3.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateCurve.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSummary.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSummary2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/calculateIC50.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateAbsoluteIC50.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/computeSensitivity.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/createIC50Summary.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/imputeData.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/imputeData2.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/scoreCurve.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/scoreCurve.L4.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/summarizeBlocks.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/summarizeSynergies.R') #complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/LL4Reverse.R') # complete
  # source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/L4Reverse.R') # complete
  options(scipen = 999) # suppress the scientific notation

  if (!all(c("block_id", "drug_row", "drug_col", "response", "conc_r", "conc_c",
             "conc_r_unit", "conc_c_unit","cell_line_id", "drug_row_id",
             "drug_col_id") %in%
           colnames(response)))
    stop("The input data must contain the following columns: ",
         "block_id, drug_row, drug_col, response,\n",
         "conc_r, conc_c, conc_r_unit, conc_c_unit, \n",
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