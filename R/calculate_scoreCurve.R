# CSS - scoreCurve
# New function used to score sensitivities given either a single-agent or a
# fixed conc (combination) columns.
# The function calculates the AUC of the log10-scaled dose-response curve.

# IMPORTANT: note that with L.4() calls, this value is already logged since the
# input concentrations are also logged.
# c1:log10(min conc) (this is the minimal nonzero concentration)
# c2:log10(max conc) (this is the maximal concentration)
# t: threshold (usually set to zero)
# own_log: log(1+10^(b*(c-x)))

scoreCurve <- function(b, c, d, m, c1, c2, t) {
  # y <- c + (d - c) / (1 + (e / x) ^ (-b)) # LL.4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)
  int_y <- (((((d - c) * own_log(-b, c2, m)) / ((-b) * log(10))) + c * c2) -
              ((((d - c) * own_log(-b, c1, m)) / ((-b) * log(10))) + c * c1))
  # int_y <- (((((a - d) * own_log(b, c, x2)) / (b * log(10))) + a * x2) -
  #            ((((a - d) * own_log(b, c, x1)) / (b * log(10))) + a * x1)) -
  #              (t * (x2 - x1))

  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}
