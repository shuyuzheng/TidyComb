# TidyComb
# Functions for calculating cell line's sensitivity to drugs or drug combination
# Copyright Shuyu Zheng
#
# Functions on this page:
#


# IC50
CalculateIC50 <- function(e, type, max.conc) {
  if (type == "LL.4") {
    ic50 <- e
  } else if (type == "L.4") {
    ic50 <- exp(e)
  }

  if (ic50 > max.conc) {
    ic50 = max(df$conc)
  }

  return (ic50)
}