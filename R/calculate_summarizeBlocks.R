# CSS - summarizeBlocks
# This function starts by summarizing the whole dataset, where each block is
# turned into one row.
# For each experiment, and for each of the two drugs in it
  # it first utilizes functions computeSensitivity() and scoreCurve() to
  # compute their DSS scores.
  # it then obtains their selected IC50 values from drugsIC50
  # it checks whether the selected IC50 values have been tested in the current
  # block or not if they were, it utilizes them with the same stated functions
  # to compute CSS scores. otherwise, it first utilizes imputeData() to predict
  # values at the missing IC50s afterwards, it utilizes the same stated
  # functions on the new data to compute CSS scores. At the end, obatined IC50
  # values and computed drug sensitivity scores are saved for the experiment at
  # hand.

summarizeBlocks <- function(raw, drugsIC50) {
  # prepare the table that will hold one row per experiment
  sensscores <- unique(raw[, c("block_id", "cell_line_name",
                               "drug_row", "drug_col")])
  rownames(sensscores) <- sensscores$block_id
  m <- sensscores$block_id
  len <- length(m)
  sensscores <- sensscores[, -1]
  sensscores$IC_row  <- NA
  sensscores$IC_col  <- NA
  sensscores$DSS_row <- NA
  sensscores$DSS_col <- NA
  sensscores$CSS_row <- NA
  sensscores$CSS_col <- NA

  options(show.error.messages = FALSE)
  # each block ID defines an experiment, i.e. a combination for a particular
  # cell line
  for (i in 1:len) {
    #progress(i,progress.bar = TRUE)
    # cat("\r", paste("Block ", i, " of ", nrow(sensscores)))
    cat("\r", i)
    # extract comb info for this block
    #print(i)
    tempmaindf <- raw[which(raw$'block_id' == m[i]), ]
    # add random noise to avoid zero variance
    tempmaindf$response <- tempmaindf$response +
                             stats::runif(length(tempmaindf$response), -0.01, 0.01)

    # DSS_c
    tempdf_c <- tempmaindf[which(tempmaindf$'conc_r' == 0),
                                 c("conc_c", "response")]
    tempdf_c <- tempdf_c[order(tempdf_c$'conc_c'), ]
    # rownames(tempdf_c) <- tempdf_c$'conc_c'
    # tempdf_c <- tempdf_c[,'response',drop = FALSE]

    #if response is unchanged after the first zero concentration, add tiny
    #cutoff to last point
    # if(var(tempdf_c$response[-1]) == 0) {
    #   tempdf_c$response[nrow(tempdf_c)] <- tempdf_c$response[nrow(tempdf_c)] +
    #                                          10^-10
    # }
    names(tempdf_c)[1] <- "conc"
    DSS_c <- computeSensitivity(tempdf_c)

    # DSS_r
    tempdf_r <- tempmaindf[which(tempmaindf$'conc_c'==0),
                           c("conc_r","response")]
    tempdf_r <- tempdf_r[order(tempdf_r$'conc_r'), ]
    # rownames(tempdf_r)=tempdf_r$'conc_r'
    # tempdf_r=tempdf_r[,'response',drop=F]
    #if response is unchanged after the first zero concentration, add tiny
    #cutoff to last point
    # if(var(tempdf_r$'response'[2:nrow(tempdf_r)])==0) {
    #  tempdf_r[nrow(tempdf_r), 'response'] <- tempdf_r[nrow(tempdf_r),
    #                                                   'response'] + 10^-10
    # }
    # if (var(tempdf_r$response[-1]) == 0) {
    #   tempdf_r$response[nrow(tempdf_r)] <- tempdf_r$response[nrow(tempdf_r)] +
    #                                        10^-10
    # }
    names(tempdf_r)[1] <- "conc"
    DSS_r = computeSensitivity(tempdf_r)

    # identify the drugs
    drug_r <- sensscores[i, 'drug_row']
    drug_c <- sensscores[i, 'drug_col']

    # pick up their selected IC50s
    ic_r <- as.numeric(drugsIC50[which(rownames(drugsIC50) == drug_r),
                                 'SelectedConcentration'])
    ic_c <- as.numeric(drugsIC50[which(rownames(drugsIC50) == drug_c),
                                 'SelectedConcentration'])

    # CSS_c
    # check if ic_r has been tested in this block
    if (ic_r %in% tempmaindf$'conc_r') {
      tempcf_c <- tempmaindf[which(tempmaindf$'conc_r' == ic_r),
                             c("conc_c","response")]
      tempcf_c <- tempcf_c[order(tempcf_c$'conc_c'), ]
      # rownames(tempcf_c) <- tempcf_c$'conc_c'
      # tempcf_c <- tempcf_c[, 'response', drop = FALSE]
    } else {
      tempcf_c <- imputeData(tempmaindf, "conc_c", "conc_r", ic_r)
    }

    # if (var(tempcf_c$response[-1]) == 0) {
    #   tempcf_c$response[nrow(tempcf_c)] <- tempcf_c$response[nrow(tempcf_c)] +
    #                                        10^-10
    # }
    names(tempcf_c)[1] <- "conc"
    CSS_c <- computeSensitivity(tempcf_c)

    # CSS_r
    # check if ic_c has been tested in this block
    if (ic_c %in% tempmaindf$'conc_c') {
      tempcf_r <- tempmaindf[which(tempmaindf$'conc_c' == ic_c),
                             c("conc_r","response")]
      tempcf_r <- tempcf_r[order(tempcf_r$'conc_r'), ]
      # rownames(tempcf_r) <- tempcf_r$'conc_r'
      # tempcf_r <- tempcf_r[, 'response', drop = FALSE]
    } else {
      tempcf_r <- imputeData(tempmaindf, "conc_r", "conc_c", ic_c)
    }

    # if (var(tempcf_r$response[-1]) == 0) {
    #   tempcf_r$response[nrow(tempcf_r)] <- tempcf_r$response[nrow(tempcf_r)] +
    #                                        10^-10
    # }
    names(tempcf_r)[1] <- "conc"
    CSS_r <- computeSensitivity(tempcf_r)

    #assign six measurements to sensscores
    sensscores[i, 4:9] <- c(ic_r, ic_c, DSS_r, DSS_c, CSS_r, CSS_c)

  }
  options(show.error.messages = TRUE)

  sensscores$CSS <- NA
  for (i in 1:nrow(sensscores)) {
    sensscores$CSS[i] <- mean(c(sensscores$CSS_r[i], sensscores$CSS_c[i]),
                              na.rm = TRUE)
  }
  return(sensscores)
}
