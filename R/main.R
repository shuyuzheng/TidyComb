# Function to take the uploaded response data
# response = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\upload\\response_template.csv", stringsAsFactors = F)
# pipeline(response)

pipeline = function(response){
  library(openxlsx)
  library(plyr)
  library(reshape2)
  library(drc)
  library(nleqslv)
  library(gplots)
  library(dplyr)
  library(SpatialExtremes)

  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_log.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_log2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData3.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/FittingSingleDrug2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateSynergy2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ZIP2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Loewe2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Bliss2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/HSA2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/BaselineCorrectionSD2.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/d1.fun.R') # used in loewe # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/d2.fun.R') # used in loewe # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/eq.R') # used in loewe # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateScore.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface2.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSurface3.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateCurve.R') # complete
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSummary.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/GenerateSummary2.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/calculateIC50.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateAbsoluteIC50.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/computeSensitivity.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/createIC50Summary.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/imputeData.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/imputeData2.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/scoreCurve.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/scoreCurve.L4.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/summarizeBlocks.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/summarizeSynergies.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/LL4Reverse.R')
  source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/L4Reverse.R')
  options(scipen = 999) # suppress the scientific notation

  if (!all(c("block_id", "drug_row", "drug_col",
             "response", "conc_r", "conc_c", "conc_r_unit", "conc_c_unit","cell_line_name") %in%
           colnames(response)))
    stop("The input data must contain the following columns: block_id, drug_row, drug_col, response,\n         conc_r, conc_c, conc_r_unit, conc_c_unit", "cell_line_name")

  print(noquote('Generating scores...'))
  response_with_scores = GenerateScore(response)

  print(noquote('Generating summary...'))
  summary = GenerateSummary2(response_with_scores)

  print(noquote('Generating surface...'))
  surface = GenerateSurface3(response_with_scores)

  print(noquote('Generating curve...'))
  curve = GenerateCurve(response_with_scores)

  print(noquote('Ready.'))
  return(list(response_with_scores=response_with_scores, summary = summary, surface = surface, curve = curve))
}