# TidyComb
# Functions for calculating drug synergy scores.
# Copyright Shuyu Zheng
#


ExtractSingleDrug <- function(response.mat, dim = "row") {
  if (dim == "row") {
    single.drug <- data.frame(response = response.mat[, "0"],
                              dose = as.numeric(rownames(response.mat)))
  } else if (dim == "col") {
    single.drug <- data.frame(response = response.mat["0", ],
                              dose = as.numeric(colnames(response.mat)))
  } else {
    stop("Values for 'dim' should be eighther 'row' or 'col'!")
  }
  rownames(single.drug) <- NULL
  return(single.drug)
}

#' Fitting single drug dose-response model
#'
#' Function \code{FitDoseResponse} fits dose-response model by using
#' \code{\link[drc]{drm}} function.
#'
#' Pre-fitting process:
#' 1. Change the 0 value in concentration into 10^-10 to avoide raising error
#' when taking log.
#' 2. If the variance of "response" values equal to 0, add 10^-10 to the last
#' "response" value.
#'
#' Model choice:
#' First use "L.4" model to fit the raw data. If error or waring occurs, use
#' "LL.4" model to fit log(raw data).
#'
#' @param data A data frame. It contains two columns:
#' \itemize{
#'   \item \strong{conc} The concentration of drugs added in experiment.
#'   \item \strong{response} The response of cell lines to drug with different
#'   concentrations.
#' }
#'
#' @param Emin A numeric or \code{NA}. It specifies the minimum value in the
#' fitted dose-response curve. Default setting is \code{NA}.
#'
#' @param Emax A numeric or \code{NA}. It specifies the maximum value in the
#' fitted dose-response curve. Default setting is \code{NA}.
#'
#' @return An object of class 'drc'. It contains imformation of fitted model.
#'
#' @export
#'
FitDoseResponse <- function (data, Emin = NA, Emax = NA) {

  if (!all(c("dose", "response") %in% colnames(data))) {
    stop('The input must contain columns: "dose", "respone".')
  }

  # nonzero concentrations to take the log
  data$dose[which(data$dose == 0)] <- 10^-10

  # ???
  if(nrow(data) != 1 & stats::var(data$response) == 0) {
    data$response[nrow(data)] <- data$response[nrow(data)] + 10 ^ -10
  }

  drug.model <- tryCatch({
    drc::drm(response ~ dose, data = data,
             fct = drc::LL.4(fixed = c(NA, Emin, Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, warning = function(w) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, Emin, Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, error = function(e) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, Emin, Emax, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  })
  return(drug.model)
}



CalculateZIP <- function(response.mat, correction = TRUE,
                         Emin = NA, Emax = NA) {
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                    Emin = Emin, Emax = Emax)))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                    Emin = Emin, Emax = Emax)))

  if (correction) {
    baseline <- (min(as.numeric(drug.row.fit)) +
                   min(as.numeric(drug.col.fit))) / 2
    response.mat <- response.mat - ((100 - response.mat) / 100 * baseline)
  }


  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                     Emin = Emin, Emax = Emax)))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                     Emin = Emin, Emax = Emax)))


  # updated.single.mat[1, -1] <- drug.col.fit[-1]
  # updated.single.mat[-1, 1] <- drug.row.fit[-1]
  n.row <- nrow(response.mat)
  n.col <- ncol(response.mat)

  # print('row')
  tmp <- data.frame(dose = as.numeric(rownames(response.mat))[-1])
  updated.col.mat <- matrix(nrow = n.row - 1, ncol = n.col - 1)

  for (i in 2:n.col) {
    # nonzero concentrations to take the log
    tmp$response <- response.mat[-1, i]

    # ????: | or &, + or -
    # if (stats::var(tmp$response, na.rm = TRUE) == 0 |
    #                length(tmp$response) == 1) {
    #   tmp$response[1] <- tmp$response[1] - 10^-10
    # }

    if(nrow(tmp) == 1) {
      # # no fitting
      fitted.response <- tmp$response - 10 ^ -10
    } else {
      tmp.min <- drug.col.fit[i]
      tmp.model <- FitDoseResponse(data = tmp, Emin = tmp.min, Emax = 100)
      fitted.response <- suppressWarnings(stats::fitted(tmp.model))
    }

    # if (fitted.inhibition[length(fitted.inhibition)] < 0)
    #  fitted.inhibition[length(fitted.inhibition)] <- tmp.min
    updated.col.mat[, i - 1] <- fitted.response
  }

  # print('col')
  tmp <- data.frame(dose = as.numeric(colnames(response.mat))[-1])
  updated.row.mat <- matrix(nrow = n.row - 1, ncol = n.col - 1)

  for (i in 2:n.row) {
    # nonzero concentrations to take the log
    tmp$response <- response.mat[i, -1]

    if(nrow(tmp) == 1) {
      # # no fitting
      fitted.response <- tmp$response - 10 ^ -10
    } else {
      tmp.min <- drug.row.fit[i]
      tmp.model <- FitDoseResponse(data = tmp, Emin = tmp.min, Emax = 100)
      fitted.response <- suppressWarnings(stats::fitted(tmp.model))
    }

    # if (fitted.inhibition[length(fitted.inhibition)] < 0)
    #  fitted.inhibition[length(fitted.inhibition)] <- tmp.min
    updated.row.mat[i - 1 , ] <- fitted.response
  }

  fitted.mat <- (updated.col.mat + updated.row.mat) / 2

  zip.mat <- matrix(nrow = n.row -1, ncol = n.col - 1)
  for (i in 1:(n.row - 1)) {
    for (j in 1:(n.col - 1)) {
      zip.mat[i, j] <- drug.row.fit[i + 1] + drug.col.fit[j + 1] -
        drug.row.fit[i + 1] * drug.col.fit[j + 1] / 100
    }
  }

  delta.mat <- fitted.mat - zip.mat

  # add 0 to first column and first row
  delta.mat <- rbind(rep(0, times = n.row -1), delta.mat)
  delta.mat <- cbind(rep(0, times = n.col), delta.mat)

  # add colname and rowname for delta.mat
  colnames(delta.mat) <- colnames(response.mat)
  rownames(delta.mat) <- rownames(response.mat)
  return(delta.mat)
}