# Function to generate the curve table
#

# response_with_scores must contain the following columns
# "block_id","conc_r","conc_c", "response", "drug_row", "drug_col",
# "conc_r_unit", "conc_c_unit", "cell_line_name", "row", "col", "synergy_zip",
# "synergy_hsa", "synergy_bliss", "synergy_loewe")

GenerateCurve <- function(response_with_scores) {
  m <- unique(response_with_scores$block_id)
  info <- dplyr::select(response_with_scores, block_id, drug_row, drug_col,
                        conc_r_unit, conc_c_unit) %>%
    unique()
  row.names(info) <- info$block_id

  curve <- list()

  options(show.error.messages = F)
  for(i in m){
    print(i)
    # cat('\r', i)
    # single drug curve fitting
    mat3 <- reshape2::acast(response_with_scores, conc_r ~ conc_c,
                            value.var = "response",
                            function(x) mean(x, na.rm = T),
                            subset = plyr::.(block_id == i))
    # mat4 <- stats::setNames(reshape2::melt(mat3),
    #                         c('conc_r', 'conc_c', 'response'))

    # mat5 <- mat4[which(mat4$conc_r==0 | mat4$conc_c==0),]
    # single.row <- mat5[which(mat5$conc_c==0),]
    single.row <- cbind(as.numeric(rownames(mat3)), mat3[, "0"])

    # nonzero concentrations to take the log
    single.row[which(single.row[, 1] == 0), 1]<- 10^-10

    # single.col <- mat5[which(mat5$conc_r==0),]
    single.col <- cbind(as.numeric(colnames(mat3)), mat3["0",])
    # nonzero concentrations to take the log
    # single.col$conc_c[which(single.col$conc_c == 0)] <- 10^-10
    single.col[which(single.col[, 1] == 0), 1]<- 10^-10

    coef.row <- tryCatch ({
      tmp <- drc::drm(single.row[, 2] ~ single.row[, 1],
                      fct = drc::LL.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
      # y1 = (tmp$coefficients[3] + tmp$coefficients[2] *
      #      (x / tmp$coefficients[4]) ^ tmp$coefficients[1]) /
      #      (1 + (x/tmp$coefficients[4]) ^ tmp$coefficients[1])
      # y2 = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2]) /
      #     (1 + exp(tmp$coefficients[1] * (log(x) - log(tmp$coefficients[4]))))
    }, warning = function(cond) {
      tmp <- drc::drm(single.row[, 2] ~ log(single.row[, 1]),
                      fct = drc::L.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, error = function(cond) {
      tmp <- drc::drm(single.row[, 2] ~ log(single.row[, 1]),
                      fct = drc::L.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
      # y1 = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2]) /
      #      (1 + exp(tmp$coefficients[1] * (x - tmp$coefficients[4])))
    }
  )

  coef.col = tryCatch ({
    tmp = drc::drm(single.col[, 2] ~ single.col[, 1], fct = drc::LL.4(),
                   control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, warning = function(cond){
      tmp = drc::drm(single.col[, 2] ~ log(single.col[, 1]),
                     fct = drc::L.4(),
                     control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, error = function(cond){
      tmp = drc::drm(single.col[, 2] ~ log(single.col[, 1]),
                     fct = drc::L.4(),
                     control = drc::drmc(errorm = FALSE, noMessage = TRUE))
      # y = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2]) /
      #     (1 + exp(tmp$coefficients[1] * (log(x) - log(tmp$coefficients[4]))))
    }
  )


  mat6 <- data.frame(matrix(NA, nrow = 2, ncol = 10))
  mat6 <- stats::setNames(mat6, c("block_id","drug_row", "drug_col",
                                  "conc_r_unit","conc_c_unit", "b", "c", "d",
                                  "e", "model"))
#
#   mat6[1,] <- c(data$drug.pairs$blockIDs[index],
#                 data$drug.pairs$drug.row[index], NA,
#                 data$drug.pairs$concRUnit[index], NA,
#                 coef.row$coefficients, as.character(coef.row$call$fct))
#   mat6[2,] <- c(data$drug.pairs$blockIDs[index], NA,
#                 data$drug.pairs$drug.col[index], NA,
#                 data$drug.pairs$concCUnit[index],
#                 coef.col$coefficients, as.character(coef.col$call$fct))
  mat6[1, ] <- c(i, info[as.character(i), "drug_row"], NA,
                 info[as.character(i), "conc_r_unit"], NA,
                 coef.row$coefficients, as.character(coef.row$call$fct))
  mat6[2, ] <- c(i, NA, info[as.character(i), "drug_col"],
                 NA, info[as.character(i), "conc_c_unit"],
                 coef.row$coefficients, as.character(coef.row$call$fct))
  curve[[i]] <- mat6

  }
  options(show.error.messages = TRUE)
  curve_table <- do.call(rbind, curve)
  return(curve_table)
}