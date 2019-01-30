#' Reshape data into response matrix
#'
#' Function \code{ReshapeBlock} transforms response data in form of molten data
#' frame into "response matrix" format.
#'
#' \code{ReshapeBlock} extract response data of one combination block and tidy
#' them into a matrix whose columns and rows refer to consentrations of two
#' drugs and the cells were filled with \em{response value} referring to
#' inhibition rate of combination.
#'
#' User also need to declear the type of "response value" in input data which is
#' defaultly set as \code{inhibition}. If the \code{type} is \code{viability},
#' \code{ReshapeBlock} will automatically transform them into \em{inhibition
#' value} by substracting them from 100.
#'
#' @param block A data frame which contains the drug combination response data
#' and it must contain the following columns: block_id, response, conc_r,
#' conc_c, conc_r_unit, conc_c_unit, cell_line_name, drug_row, drug_col.
#'
#' @return A list contains "response matrix" for each pair of drugs.
#'
#' @export
#'
ReshapeBlock <- function(block) {

  response.mat <- reshape2::acast(block, conc_r ~ conc_c,
                                  value.var = "response",
                                  function(x) mean(x, na.rm = TRUE))


  return(response.mat)
}

ReshapeData2 <- function(data, type = "inhibition") {

  id.drug.comb <- unique(data$block_id)
  dose.response.mats <- list()
  drug.pairs <- data.frame(drug.row = character(length(id.drug.comb)),
                           drug.col = character(length(id.drug.comb)),
                           concRUnit = character(length(id.drug.comb)),
                           concCUnit = character(length(id.drug.comb)),
                           blockIDs = numeric(length(id.drug.comb)),
                           stringsAsFactors = FALSE)
  for (i in 1:length(id.drug.comb)) {
    tmp.mat <- data[which(data$block_id == id.drug.comb[i]),]
    if (type == "viability") {
      tmp.mat$Inhibition <- 100 - tmp.mat$response
    } else {
      tmp.mat$Inhibition <- tmp.mat$response
    }
    conc.col <- tmp.mat$conc_c[which(tmp.mat$row == 1)]
    conc.col <- conc.col[order(tmp.mat$col[which(tmp.mat$row == 1)])]
    conc.row <- tmp.mat$conc_r[which(tmp.mat$col == 1)]
    conc.row <- conc.row[order(tmp.mat$row[which(tmp.mat$col == 1)])]
    response.mat <- reshape2::acast(tmp.mat, conc_r ~ conc_c,
                                    value.var = "Inhibition",
                                    function(x) mean(x, na.rm = TRUE))
    # replicates are averaged

    if (which.max(conc.row) == 1 & which.max(conc.col) == 1) {
      response.mat <- t(apply(apply(response.mat, 2, rev), 1, rev))
    } else if (which.max(conc.row) == length(conc.row) &
               which.max(conc.col) == 1) {
      response.mat <- t(apply(response.mat, 1, rev))
    } else if (which.max(conc.row) == 1 &
               which.max(conc.col) == length(conc.col)) {
      response.mat <- apply(response.mat, 2, rev)
    }
    conc.runit <- unique(tmp.mat$conc_r_unit)
    conc.cunit <- unique(tmp.mat$conc_c_unit)
    drug.row <- unique(tmp.mat$drug_row)
    drug.col <- unique(tmp.mat$drug_col)
    drug.pairs$drug.row[i] <- drug.row
    drug.pairs$drug.col[i] <- drug.col
    drug.pairs$concRUnit[i] <- conc.runit
    drug.pairs$concCUnit[i] <- conc.cunit
    dose.response.mats[[i]] <- response.mat
  }
  drug.pairs$blockIDs <- id.drug.comb
  return(list(dose.response.mats = dose.response.mats, drug.pairs = drug.pairs))
}
