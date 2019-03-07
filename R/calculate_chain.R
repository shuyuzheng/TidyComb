# TidyComb
# Functions for pipe all calculation togather
# Copyright Shuyu Zheng
#
# Functions on this page:
#

MergeSynergy <- function(response.mat){
  zip <- reshape2:melt(CalculateZIP(response.mat),
                       valuename = "synergy_zip")
  bliss <- reshape2:melt(CalculateBliss(response.mat),
                         valuename = "synergy_bliss")
  hsa <- reshape2:melt(CalculateHSA(response.mat),
                       valuename = "synergy_hsa")
  loewe <- reshape2:melt(CalculateLoewe(response.mat),
                         valuename = "synergy_loewe")
  synergy <- merge(zip)
}