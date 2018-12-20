# CSS - createIC50Summary
# This function takes in the raw data as input
# It goes over every block
# For each of the block's 2 drugs, a drm function is fit and the obtained IC50
# extracted using function calculateIC50(). The IC50 is saved in the drugsIC50
# table's first column in the row corresponding to the drug at hand.
# The drugsIC50 table's second column will hold one IC50 per drug (being the
# mean of all obtained IC50s for that drug)
#
# IMPORTANT NOTE: the calculation of one IC50 per drug is independent of the
# specific cell lines from which the vector of IC50 values where obtained for
# that drug. This will permit the comparability of CSS scores later across cell
# lines (rather than having them as merely relative measures). Nevertheless, the
# current calculation is still dataset-dependent, rendering CSS values
# comparable ONLY WITHIN THE SAME DATASET. After a mean IC50 is defined per
# drug, the set of all tested concentrations for that drug in that dataset is
# collected and stored. A selected concentration per drug is finally determined,
# as being the closest tested concentration to the mean IC50. A summary table
# of the above is returned, with the first column dropped

createIC50Summary <- function(raw) {
  # blocks <- sort(unique(raw$block_id)) # sort added
  blocks <- unique(raw$block_id)
  len <- length(blocks)
  drugs <- unique(c(raw$drug_row, raw$drug_col))

  # create a data frame that will hold IC50 values
  drugsIC50 <- data.frame(matrix(NA, ncol = 4, nrow = length(drugs)))
  rownames(drugsIC50) <- drugs
  colnames(drugsIC50) <- c("IC50values", "MeanIC50", "TestedConcentrations",
                           "SelectedConcentration")

  # go over each block (i.e. each combination)
  options(show.error.messages = FALSE)
  for (i in 1:len) {
    #progress(i,progress.bar = TRUE)
    # cat("\r", paste("Block ", i, " of ", length(blocks)))
    cat("\r", i)
    tempdf <- raw[which(raw$block_id == blocks[i]), ]
    # tempd_r <- unique(tempdf$drug_row)
    # tempd_c <- unique(tempdf$drug_col)

    tempd_r <- tempdf$drug_row[1] # drug row
    tempd_c <- tempdf$drug_col[1] # drug col
    # tempc_r <- paste(sort(unique(tempdf$conc_r)), collapse = ";")
    # tempc_c <- paste(sort(unique(tempdf$conc_c)), collapse = ";")
    #print (paste(i, tempd_r, tempc_r, tempd_c, tempc_c), sep=" ")

    #for tempd_c
    tempdf_c <- tempdf[which(tempdf$conc_r==0), c("conc_c", "response")]
    tempdf_c <- tempdf_c[order(tempdf_c$conc_c), ]
    # rownames(tempdf_c) <- tempdf_c$conc_c
    # tempdf_c <- tempdf_c[,"response",drop = FALSE]

    # Call calculateIC50()
    names(tempdf_c)[1] <- "conc"
    # ic_c <- calculateIC50(tempdf_c)
    ic_c <- CalculateAbsoluteIC50(tempdf_c)

    # if (is.na(drugsIC50[tempd_c, "IC50values"])) {
    #   drugsIC50[tempd_c, "IC50values"] <- ic_c
    # } else {
    #   drugsIC50[tempd_c, "IC50values"] <- paste(drugsIC50[tempd_c,
    #                                                       "IC50values"],
    #                                             ic_c, sep = ";")
    # }
    drugsIC50[tempd_c, "IC50values"] <- paste(drugsIC50[tempd_c, "IC50values"],
                                              ic_c, sep = ";")

    #for tempd_r
    tempdf_r <- tempdf[which(tempdf$conc_c == 0), c("conc_r", "response")]
    tempdf_r <- tempdf_r[order(tempdf_r$conc_r), ]
    # rownames(tempdf_r) <- tempdf_r$conc_r
    # tempdf_r <- tempdf_r[, "response", drop = FALSE]

    # Call calculateIC50()
    names(tempdf_r)[1] <- "conc"
    ic_r <- CalculateAbsoluteIC50(tempdf_r)
    # if (is.na(drugsIC50[tempd_r, "IC50values"])){
    #   drugsIC50[tempd_r, "IC50values"] <- ic_r
    # } else {
    #   drugsIC50[tempd_r, "IC50values"] <- paste(drugsIC50[tempd_r,
    #                                                       "IC50values"],
    #                                             ic_r, sep = ";")
    # }
    drugsIC50[tempd_r, "IC50values"] <- paste(drugsIC50[tempd_r, "IC50values"],
                                              ic_r, sep = ";")
  }
  options(show.error.messages = TRUE)

  # calculate the average IC50 per drug and store it in column 2
  for (i in 1:nrow(drugsIC50)) {
    tempv <- as.numeric(unlist(strsplit(drugsIC50$"IC50values"[i],
                                        split = ";")))
    # drugsIC50$"MeanIC50"[i] <- mean(tempv, na.rm=TRUE)
    drugsIC50$"MeanIC50"[i] <- stats::median(tempv, na.rm=TRUE)
  }

  # obtain the list of concentrations tested per drug
  # ! eliminate the zero concentration from the set so that it is never selected as a drug's IC50
  for (i in drugs) {
    concset <- union(raw$conc_r[which(raw$drug_row == i)],
                     raw$conc_c[which(raw$drug_col == i)])
    concset <- unique(concset)
    concset <- sort(concset)
    if (as.numeric(concset[1]) == 0) {
      concset <- concset[-1]
    }
    drugsIC50[which(rownames(drugsIC50) == i),
              "TestedConcentrations"] <- paste(concset, collapse = "; ")
  }

  # # define the selected concentration as the one closest to meanIC50
  # for (i in 1:nrow(drugsIC50)) {
  #   concset <- as.numeric(unlist(strsplit(drugsIC50$TestedConcentrations[i],
  #                         split = "; ")))
  #   drugsIC50$SelectedConcentration[i] <- concset[which.min(abs(concset -
  #                                         as.numeric(drugsIC50$MeanIC50[i])))]
  # }

  drugsIC50$SelectedConcentration <- drugsIC50$MeanIC50

  # # drop useless first column
  # drugsIC50 <- drugsIC50[, -1]
  return(drugsIC50)
}
