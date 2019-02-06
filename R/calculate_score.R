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
#' @param ...
#'
#' @return An object of class 'drc'. It contains imformation of fitted model.
#'
#' @export
#'
FitDoseResponse <- function (data, ...) {

  if (!all(c("dose", "response") %in% colnames(data))) {
    stop('The input must contain columns: "dose", "respone".')
  }

  # nonzero concentrations to take the log
  data$dose[which(data$dose == 0)] <- 10^-10

  # ???
  if(nrow(data)!= 1 & stats::var(data$response) == 0) {
    data$response[nrow(data)] <- data$response[nrow(data)]
                                            + 10^-10
  }

  drug.model <- tryCatch({
    drc::drm(response ~ dose, data = data,
             fct = drc::LL.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, warning = function(w) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  }, error = function(e) {
    drc::drm(response ~ log(dose), data = data,
             fct = drc::L.4(fixed = c(NA, NA, NA, NA)),
             na.action = stats::na.omit,
             control = drc::drmc(errorm = FALSE, noMessage = TRUE))
  })
  return(drug.model)
}



CalculateZIP <- function(response.mat, correction = TRUE,
                         Emin = NA, Emax = NA) {
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                fixed = c(NA, Emin, Emax, NA))))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                fixed = c(NA, Emin, Emax, NA))))

  if (correction) {
    baseline <- (min(as.numeric(drug.row.fit)) +
                   min(as.numeric(drug.col.fit))) / 2
    response.mat <- response.mat - ((100 - response.mat) / 100 * baseline)
  }


  # updated.single.mat[1, -1] <- drug.col.fit[-1]
  # updated.single.mat[-1, 1] <- drug.row.fit[-1]
  n.row <- nrow(response.mat)
  n.col <- ncol(response.mat)

  # print('row')
  tem <- data.frame(dose = as.numeric(rownames(response.mat))[-1])
  updated.col.mat <- matrix(nrow = n.row - 1, ncol = n.col - 1)




  if (nrow())
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
      fitted.inhibition <- tmp$inhibition
    } else {
      tmp.min <- drug.col.fit[i]
      tmp.model <- FitDoseResponse(data = tmp, fixed = c(NA, tmp.min, 100, NA))
      fitted.inhibition <- suppressWarnings(stats::fitted(tmp.model))
    }

    # if (fitted.inhibition[length(fitted.inhibition)] < 0)
    #  fitted.inhibition[length(fitted.inhibition)] <- tmp.min
    updated.col.mat[, i] <- fitted.inhibition
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
      fitted.inhibition <- tmp$inhibition
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