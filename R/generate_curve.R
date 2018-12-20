# Function to generate the curve table
#

# response_with_scores must contain the following columns
# c("block_id","conc_r","conc_c", "response", "drug_row", "drug_col", "conc_r_unit",
#      "conc_c_unit", "cell_line_name", "row", "col", "synergy_zip", "synergy_hsa", "synergy_bliss", "synergy_loewe")

GenerateCurve <- function(response_with_scores) {
  m <- unique(response_with_scores$block_id)
  len <- length(m)
  data <- ReshapeData2(response_with_scores, data.type = "inhibition")
  curve <- list()

  options(show.error.messages = F)
  for(i in 1:len){
    print(m[i])
    # cat('\r', i)
    # single drug curve fitting
    index <- which(data$drug.pairs$blockIDs == m[i])
    mat3 <- data$dose.response.mats[[index]]
    mat4 <- stats::setNames(reshape2::melt(mat3), c('conc_r', 'conc_c', 'response'))

    mat5 <- mat4[which(mat4$conc_r==0 | mat4$conc_c==0),]
    single.row <- mat5[which(mat5$conc_c==0),]
    # nonzero concentrations to take the log
    single.row$conc_r[which(single.row$conc_r == 0)] <- 10^-10
    single.col <- mat5[which(mat5$conc_r==0),]
    # nonzero concentrations to take the log
    single.col$conc_c[which(single.col$conc_c == 0)] <- 10^-10

    coef.row <- tryCatch ({
      tmp <- drc::drm(single.row$response ~ single.row$conc_r,
                      fct = drc::LL.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
      # y1 = (tmp$coefficients[3] + tmp$coefficients[2] *
      #      (x / tmp$coefficients[4]) ^ tmp$coefficients[1]) /
      #      (1 + (x/tmp$coefficients[4]) ^ tmp$coefficients[1])
      # y2 = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2]) /
      #     (1 + exp(tmp$coefficients[1] * (log(x) - log(tmp$coefficients[4]))))
    }, warning = function(cond) {
      tmp <- drc::drm(single.row$response ~ log(single.row$conc_r),
                      fct = drc::L.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, error = function(cond) {
      tmp <- drc::drm(single.row$response ~ log(single.row$conc_r),
                      fct = drc::L.4(),
                      control = drc::drmc(errorm = FALSE, noMessage = TRUE))
      # y1 = tmp$coefficients[2] + (tmp$coefficients[3] - tmp$coefficients[2]) /
      #      (1 + exp(tmp$coefficients[1] * (x - tmp$coefficients[4])))
    }
  )

  coef.col = tryCatch ({
    tmp = drc::drm(single.col$response ~ single.col$conc_c, fct = drc::LL.4(),
                   control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, warning = function(cond){
      tmp = drc::drm(single.col$response ~ log(single.col$conc_c),
                     fct = drc::L.4(),
                     control = drc::drmc(errorm = FALSE, noMessage = TRUE))
    }, error = function(cond){
      tmp = drc::drm(single.col$response ~ log(single.col$conc_c),
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

  mat6[1,] <- c(data$drug.pairs$blockIDs[index],
                data$drug.pairs$drug.row[index], NA,
                data$drug.pairs$concRUnit[index], NA,
                coef.row$coefficients, as.character(coef.row$call$fct))
  mat6[2,] <- c(data$drug.pairs$blockIDs[index], NA,
                data$drug.pairs$drug.col[index], NA,
                data$drug.pairs$concCUnit[index],
                coef.col$coefficients, as.character(coef.col$call$fct))
  curve[[i]] <- mat6

  }
  options(show.error.messages = TRUE)
  curve_table <- do.call(rbind, curve)
  return(curve_table)
}