# CSS - scoreCurve for L.4() model
# New function used to score sensitivities given either a single-agent or a
# fixed conc (combination) columns.
# The function calculates the AUC of the log10-scaled dose-response curve
# own_log2: log(1+exp(x))

scoreCurve.L4 <- function(b, c, d, e, c1, c2, t) {
  # y <- c + (d - c) / (1 + exp(b * (x - e))) # L4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)

  int_y <- d * (c2 - c1) + ((c - d) / b) *
           (own_log2(b * (c2 - e)) - own_log2(b * (c1 - e)))
  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}
