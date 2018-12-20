# CSS - summarizeSynergies
# This function receives a data frame containing synergy values per experiment
# Per experiment, there are as many synergy values as their are
# drug-concentrations' combination pairs (excluding zeroes)
# times four different types of synergies
# The function simply equates a mean per synergy type and per experiment, so
# that one block is represented by one row

summarizeSynergies <- function(syn) {
  synscores <- data.frame(matrix(NA, ncol = 4,
                                 nrow = length(unique(syn$'block_id'))))
  # rownames(synscores) <- unique(syn$'block_id')
  colnames(synscores) <- c("synergy_zip", "synergy_hsa",
                           "synergy_bliss", "synergy_loewe")
  m <- unique(syn$'block_id')
  len <- length(m)
  for (i in 1:len) {
    # progress(i,progress.bar = TRUE)
    # cat("\r",paste("Block ", i, " of ", length(unique(syn$'block_id'))))
    # cat("\r", i)
    tempmat <- syn[which(syn$'block_id' == m[i]), ]
    synscores[i, ] <- c(mean(tempmat$'synergy_zip', na.rm = TRUE),
                        mean(tempmat$'synergy_hsa',na.rm = TRUE),
                        mean(tempmat$'synergy_bliss',na.rm = TRUE),
                        mean(tempmat$'synergy_loewe',na.rm = TRUE))
  }
  rownames(synscores) <- m
  return(synscores)
}
