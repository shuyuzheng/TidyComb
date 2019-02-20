# TidyComb
# Functions for calculating drug synergy scores.
# Copyright Shuyu Zheng
#


#' Extract single drug
#'
#' \code{ExtractSingleDrug} extracts the dose-response values of single drug (
#' drug added in column or row) from a drug combination dose-response matrix.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @param dim A character. It should be either "col" or "row" to indicate which
#' drug's dose-response value will be extracted.
#'
#' @return A data frame. It contains two variables:
#' \itemize{
#' \item \strong{dose} The concertration of drug.
#' \item \strong{response} The cell's response (inhibation rate) to
#' corresponding drug concertration.
#' }
#'
#' @export
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

#' Base line correction
#'
#' \code{CorrectBaseLine} adjusts the base line of drug combination
#' dose-response matrix up to positive values.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @return A matrix which base line have been adjusted.
#' @export
CorrectBaseLine <- function(response.mat){
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                  Emin = Emin, Emax = Emax)))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                  Emin = Emin, Emax = Emax)))

  baseline <- (min(as.numeric(drug.row.fit)) +
                 min(as.numeric(drug.col.fit))) / 2

  response.mat <- response.mat - ((100 - response.mat) / 100 * baseline)
  return(response.mat)
}

#' Calculate ZIP synergy score
#'
#' \code{CalculateZIP} calculates the \eqn{\Delta} score matrix for a block of
#' drug combination by using Zero Interaction Potency (ZIP) method.
#'
#' zero interaction potency (ZIP) is a reference model for evaluating the
#' interaction between two drugs. It captures the drug interaction relationships
#' by comparing the change in the potency of the dose-response curves between
#' individual drugs and their combinations. More details about this model could
#' be found in original publication
#' \href{10.1016/j.csbj.2015.09.001}{(Yadav.et.al., 2015)}.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @param correction A logical value. It indicates whether \emph{baseline
#' correction} needed before calculation. Default value is \code{TRUE}, which
#' means the correction will be done.
#'
#' @param Emin A numeric or \code{NA}. It specifies the minimum value in the
#' fitted dose-response curve. Default setting is \code{NA}.
#'
#' @param Emax A numeric or \code{NA}. It specifies the maximum value in the
#' fitted dose-response curve. Default setting is \code{NA}.
#'
#' @return A matrix with  \eqn{\Delta} score calculated via Zero Interaction
#' Potency (ZIP) method
#'
#' @export
CalculateZIP <- function(response.mat, # correction = TRUE,
                         Emin = NA, Emax = NA) {
  # if (correction) {
  #   response.mat <- CorrectBaseLine(response.mat)
  # }

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
    drug.row <- ExtractSingleDrug(response.mat, dim = "row")
    drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                    Emin = Emin, Emax = Emax)))

    drug.col <- ExtractSingleDrug(response.mat, dim = "col")
    drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                    Emin = Emin, Emax = Emax)))
  }

  n.row <- nrow(response.mat)
  n.col <- ncol(response.mat)

  # generate drug_row fitting matrix
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

  # generate drug_col fitting matrix
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

  # clean up
  gc()
}

#' Calculate Bliss synergy score
#'
#'
#'
#' @param response.mat
#' @param correction
#' @param Emin
#' @param Emax
#'
#' @return
#' @export
#'
#' @examples
CalculateBliss <- function (response.mat, # correction = TRUE,
                            Emin = NA, Emax = NA) {
  # if (correction) {
  #   response.mat <- CorrectBaseLine(response.mat)
  # }
  #
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- drug1.response[i] + drug2.response[j] -
        drug1.response[i] * drug2.response[j]/100
    }
  }
  syn.mat <- response.mat - ref.mat
  return(syn.mat)
}