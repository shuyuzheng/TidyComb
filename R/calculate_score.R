# TidyComb
# Functions for calculating drug synergy scores.
# Copyright Shuyu Zheng
#

#' Fitting single drug dose-response model
#'
#' Function \code{FittingSingleDrug} fits dose-response model by using
#' \code{\link[drc]drm} function.
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
#' @param response A matrix. It contains two columns:
#' \itemize{
#'   \item \strong{conc} The concentration of drugs added in experiment.
#'   \item \strong{response} The response of cell lines to drug with different
#'   concentrations.
#' }
#'
#' @param ...
#'
#' @return A list. It contains two elements:
#' \itemize{
#'   \item \strong{fitted} Fitted values.
#'   \item \strong{model} Fitted model.
#' }
#'
#' @export
#'
FittingSingleDrug <- function (response, ...) {

  if (!all(c("conc", "response") %in% colnames(response))) {
    stop('The input must contain columns: "conc", "respone".')
  }

  # nonzero concentrations to take the log
  response[which(response[, "conc"] == 0), "conc"] <- 10^-10

  # ???
  if(nrow(response)!= 1 & stats::var(response[, "response"]) == 0) {
    response[nrow(response), "response"] <- response[nrow(response), "response"]
                                            + 10^-10
  }

  drug.model <- tryCatch({
    drc::drm(response ~ conc, data = response,
             fct = drc::LL.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, warning = function(w) {
    drc::drm(response ~ log(conc), data = response,
             fct = drc::L.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, error = function(e) {
    drc::drm(response ~ log(conc), data = response,
             fct = drc::L.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  })
  drug.fitted <- suppressWarnings(stats::fitted(drug.model))

  res <- list(fitted() = drug.fitted,
              model = drug.model)
  return(res)
}

ExtractSingleDrug <- function(response.mat, dim = "row") {
  if (dim == "row") {
    single.drug <- cbind(rownames(response.mat), response.mat[, "0"])
    colnames(single.drug) <- c("conc", "response")
  } else if (dim = "col") {
    single.drug <- cbind(colnames(response.mat), response.mat["0", ])
  } else {
    stop("Values for 'dim' should be eighther 'row' or 'col'!")
  }
}

#' Correct baseline of fitted drug dose-response curve
#'
#' @param response A data frame. It contains the response data for one pair of
#' drugs. It must contain colums
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
BaselineCorrectionSD <- function (response.mat, drug.row.fitted,
                                  drug.col.fitted, ...) {

  return(corrected.mat)
}

CalculateZIP <- function(response.mat, correction = TRUE,
                         Emin = NA, Emax = NA) {
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fitted <- FittingSingleDrug(drug.row, fixed = c(NA, Emin, Emax, NA))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fitted <- FittingSingleDrug(drug.col, fixed = c(NA, Emin, Emax, NA))

  if (correction) {
    baseline <- (min(as.numeric(drug.row.fitted$fitted)) +
                   min(as.numeric(drug.col.fitted$fitted)))/2
    response.mat <- response.mat - ((100 - response.mat)/100 * baseline)
  }

  drug.col.response <- single.fitted$drug.col.fitted # first row
  drug.row.response <- single.fitted$drug.row.fitted # first column

  updated.single.mat <- mat.or.vec(nrow(response.mat), ncol(response.mat))
  colnames(updated.single.mat) <- colnames(response.mat)
  rownames(updated.single.mat) <- rownames(response.mat)

  updated.single.mat[1, c(2:ncol(response.mat))] <- drug.col.response[-1]
  updated.single.mat[c(2:nrow(response.mat)), 1] <- drug.row.response[-1]
  updated.col.mat <- updated.single.mat

  # print('row')
  for (i in 2:ncol(response.mat)) {
    # print(i)
    tmp <- as.data.frame(mat.or.vec(nrow(response.mat) - 1, 0))
    tmp$dose <- as.numeric(rownames(response.mat))[-1]
    # nonzero concentrations to take the log
    tmp$dose[which(tmp$dose==0)] <- 10^-10
    tmp$inhibition <- response.mat[c(2:nrow(response.mat)), i]

    tmp.min <- updated.single.mat[1, i]
    if (stats::var(tmp$inhibition, na.rm = TRUE) == 0 |
                   length(tmp$inhibition) == 1) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }

    if(nrow(tmp) == 1) {
      # # no fitting
      fitted.inhibition = response.mat[c(2:nrow(response.mat)), i]
    } else {
      tmp.model = tryCatch(
        {
          tmp.model <- drc::drm(inhibition ~ dose, data = tmp,
                           fct = drc::LL.4(fixed = c(NA, tmp.min, 100, NA)),
                           na.action = stats::na.omit,
                           control = drc::drmc(errorm = FALSE,
                                               noMessage = TRUE))
        },
        warning = function(w) {
         tmp.model <- drc::drm(inhibition ~ log(dose), data = tmp,
                          fct = drc::L.4(fixed = c(NA, tmp.min, 100, NA)),
                          na.action = stats::na.omit,
                          control = drc::drmc(errorm = FALSE,
                                              noMessage = TRUE))
        },
        error = function(e) {
          tmp.model <- drc::drm(inhibition ~ log(dose), data = tmp,
                           fct = drc::L.4(fixed = c(NA, tmp.min, 100, NA)),
                           na.action = stats::na.omit,
                           control = drc::drmc(errorm = FALSE,
                                               noMessage = TRUE))

        }
      )
      fitted.inhibition = suppressWarnings(stats::fitted(tmp.model))
    }

    tmp$fitted.inhibition <- fitted.inhibition
    # if (tmp$fitted.inhibition[nrow(response.mat) - 1] < 0)
    #  tmp$fitted.inhibition[nrow(response.mat) - 1] <- tmp.min
    updated.col.mat[c(2:nrow(response.mat)), i] <- tmp$fitted.inhibition
  }

  # print('col')
  updated.row.mat <- updated.single.mat
  for (i in 2:nrow(response.mat)) {
    # print(i)
    tmp <- as.data.frame(mat.or.vec(ncol(response.mat) - 1, 0))
    tmp$dose <- as.numeric(colnames(response.mat))[-1]
    # nonzero concentrations to take the log
    tmp$dose[which(tmp$dose==0)] =  10^-10

    tmp$inhibition <- response.mat[i, c(2:ncol(response.mat))]
    tmp.min <- updated.single.mat[i, 1]

    if (stats::var(tmp$inhibition, na.rm = TRUE) == 0 |
        length(tmp$inhibition) == 1) {
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }

    if(nrow(tmp)==1) {
      fitted.inhibition = response.mat[i, c(2:ncol(response.mat))]
    } else {
      tmp.model = tryCatch(
        {
          tmp.model <- drc::drm(inhibition ~ dose, data = tmp,
                               fct = drc::LL.4(fixed = c(NA, tmp.min, 100, NA)),
                               na.action = stats::na.omit,
                               control = drc::drmc(errorm = FALSE,
                                                   noMessage = TRUE))
        },
        warning = function(w) {
          tmp.model <- drc::drm(inhibition ~ log(dose), data = tmp,
                                fct = drc::L.4(fixed = c(NA, tmp.min, 100, NA)),
                                na.action = stats::na.omit,
                                control = drc::drmc(errorm = FALSE,
                                                    noMessage = TRUE))

        },
        error = function(w) {
          tmp.model<-drc::drm(inhibition ~ log(dose), data = tmp,
                              fct = drc::L.4(fixed = c(NA, tmp.min, 100, NA)),
                              na.action = stats::na.omit,
                              control = drc::drmc(errorm = FALSE,
                                                  noMessage = TRUE))
        }

      )
      fitted.inhibition <- suppressWarnings(stats::fitted(tmp.model))
    }

    # if (tmp$fitted.inhibition[ncol(response.mat) - 1] < 0)
    #  tmp$fitted.inhibition[ncol(response.mat) - 1] <- tmp.min
    tmp$fitted.inhibition <- fitted.inhibition
    updated.row.mat[i, c(2:ncol(response.mat))] <- tmp$fitted.inhibition
  }

  fitted.mat <- (updated.col.mat + updated.row.mat)/2
  zip.mat <- updated.single.mat
  for (i in 2:nrow(updated.single.mat)) {
    for (j in 2:ncol((updated.single.mat))) {
      zip.mat[i, j] <- updated.single.mat[i, 1] + updated.single.mat[1, j] -
                       updated.single.mat[i, 1] * updated.single.mat[1, j]/100
    }
  }

  # fitted.mat[1, 1] <- 0
  # zip.mat[1, 1] <- 0
  # fitted.mat <- apply(fitted.mat, c(1, 2), function(x) ifelse(x > 100, 100, x))

  delta.mat <- (fitted.mat - zip.mat)
  return(delta.mat)
}