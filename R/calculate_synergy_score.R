# TidyComb
# Functions for calculating drug synergy scores.
# Copyright Shuyu Zheng
#
# Functions in this page:
#
# FitDoseResponse: Fitting single drug dose-response model
# CorrectBaseLine: Do base line correction to matrix
# CalculateZIP/Bliss/HSA/Loewe: Calculat synergy scores
# eq.LL4/L4.LL4/L4: Four functions to calculate loewe score in CalculateLoewe.
# fun: Function used in CalculateLoewe

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
CorrectBaseLine <- function(response.mat, ...){
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                                 ...)))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                                 ...)))

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
CalculateZIP <- function(response.mat, ...) {
  # if (correction) {
  #   response.mat <- CorrectBaseLine(response.mat)
  # }

  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.row,
                                                                 ...)))

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.fit <- suppressWarnings(stats::fitted(FitDoseResponse(drug.col,
                                                                 ...)))

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
#' \code{CalculateBliss} calculates the synergy score matrix for a block of
#' drug combination by using a druginteraction reference model introduced by
#' C. I. Bliss in 1939.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The expected effect of two
#' drugs acting independently". More details about this model could be found in
#' original publication:
#' \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1744-7348.1939.tb06990.x}{(Bliss, 1939)}.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @return A matrix with synergy score calculated via reference model introduced
#' by C. I. Bliss.
#'
#' @export
CalculateBliss <- function (response.mat) {
  drug.row <- response.mat[, 1]
  drug.col <- response.mat[1, ]
  reference.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      reference.mat[i, j] <- drug.row[i] + drug.col[j] -
        drug.row[i] * drug.col[j]/100
    }
  }
  synergy.mat <- response.mat - reference.mat
  return(synergy.mat)

  # clean up
  gc()
}

#' Calculate HSA synergy score
#'
#' \code{CalculateHSA} calculates the synergy score matrix for a block of
#' drug combination by using Highest Single Agent (HSA) reference model.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The reference effect of drug
#' combination is the maximal single drug effect". More details about this model
#' could be found in original publication:
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/2692037}{(Berenbaum, 1989)}.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @return A matrix with synergy score calculated via Highest Single Agent (HSA)
#' reference model.
#'
#' @export
CalculateHSA <- function(response.mat) {
  drug.row <- response.mat[, 1]
  drug.col <- response.mat[1, ]
  reference.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      reference.mat[i, j] <- max(drug.row[i], drug.col[j])
    }
  }
  synergy.mat <- response.mat - reference.mat
  return(synergy.mat)

  #clean up
  gc()
}

# Four functions to calculate loewe
eq.LL4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / (drug.col.par[4] * (((x - drug.col.par[3]) /
                            (drug.col.par[2] - x)) ^ (1/drug.col.par[1]))) +
    x2 / (drug.row.par[4] * (((x - drug.row.par[3]) /
                              (drug.row.par[2] - x)) ^ (1/drug.row.par[1]))) - 1
}# Eq.8 in the ZIP paper

eq.L4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / exp((drug.col.par[4] + log((drug.col.par[3] - x) /
                                    (x - drug.col.par[2])) / drug.col.par[1])) +
    x2 / exp((drug.row.par[4] + log((drug.row.par[3] - x) /
                                  (x - drug.row.par[2])) / drug.row.par[1])) -1
}# x1, x2 to be log scaled

eq.LL4.L4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / (drug.col.par[4] * (((x - drug.col.par[3]) /
                              (drug.col.par[2] - x)) ^ (1 / drug.col.par[1]))) +
    x2 / exp((drug.row.par[4] + log((drug.row.par[3] - x) /
                                  (x - drug.row.par[2])) / drug.row.par[1])) -1
}# x2 to be log-scaled

eq.L4.LL4 <- function(x, x1, x2, drug.col.par, drug.row.par) {
  x1 / exp((drug.col.par[4] + log((drug.col.par[3] - x) /
                                    (x - drug.col.par[2])) / drug.col.par[1])) +
    x2 / (drug.row.par[4] * (((x - drug.row.par[3]) /
                            (drug.row.par[2] - x)) ^ (1 / drug.row.par[1]))) - 1
}# x1 to be log-scaled

# function used to calculate loewe if termination code from 'nleqslv' function
# is -10, 1, or 2, which mean:
# * -10 User supplied Jacobian is most likely incorrect.
# * 1 Function criterion is near zero. Convergence of function values has been
# achieved.
# * 2 x-values within tolerance. This means that the relative distance between two
# consecutive x-values is smaller than xtol but that the function value
# criterion is still larger than ftol. Function values may not be near zero;
# therefore the user must check if function values are acceptably small.
#
fun <- function(col_conc, row_conc, drug.par, model) {
  # LL.4, conc must be raw
  if(model == "LL.4") {
    conc = col_conc + row_conc
    (drug.par[3] + drug.par[2] *
        (conc / drug.par[4]) ^ drug.par[1]) /
      (1 + (conc / drug.par[4]) ^ drug.par[1])
  } else if (model == "L.4"){# L.4, conc must be logscaled, ie. log(conc)
    conc = log(col_conc+row_conc)
    (drug.par[2] + (drug.par[3] - drug.par[2]) /
        (1 + exp(drug.par[1] * (conc - drug.par[4]))))
  } else {
    stop("Model type is incorrect. Available values are 'LL.4' or 'L.4' ")
  }
}

#' Calculate Loewe synergy score
#'
#' \code{CalculateLoewe} calculates the synergy score matrix for a block of
#' drug combination by using a druginteraction reference model introduced by
#' Loewe in 1953.
#'
#' This model is a reference model for evaluating the interaction between two
#' drugs. The basic assumption of this model is "The referece effect of drug
#' combination is the expected effect of a drug combined with itself". More
#' details about this model could be found in original publication:
#' \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1744-7348.1939.tb06990.x}{(Bliss, 1939)}.
#'
#' @param response.mat A drug cobination dose-response matrix. It's column name
#' and row name are representing the concerntrations of drug added to column and
#' row, respectively. The values in matrix indicate the inhibition rate to cell
#' growth.
#'
#' @param quiet A logical value. If it is \code{TRUE} then the warning message
#' will not show during calculation.
#'
#' @return A matrix with synergy score calculated via reference model introduced
#' by C. I. Bliss.
#'
#' @expor
CalculateLoewe <- function (response.mat, quiet = TRUE, ...) {
  if (quiet) {
    options(warn = -1)
  }
  drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  drug.row.model <- FitDoseResponse(drug.row)#, ...)
  drug.row.par <- stats::coef(drug.row.model)
  drug.row.fct <- drug.row.model$call$fct[[1]][[3]]

  drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  drug.col.model <- FitDoseResponse(drug.col)#, ...)
  drug.col.par <- stats::coef(drug.col.model)
  drug.col.fct <- drug.col.model$call$fct[[1]][[3]]

  drug.row$dose[drug.row$dose == 0] = 10^-10 # avoid log(0)
  drug.col$dose[drug.col$dose == 0] = 10^-10 # avoid log(0)

  loewe.mat <- response.mat
  eq <- switch (paste(drug.col.fct, drug.row.fct),
                "LL.4 LL.4" = eq.LL4.LL4,
                "L.4 L.4"   = eq.L4.L4,
                "LL.4 L.4"  = eq.LL4.L4,
                "L.4 LL.4"  = eq.L4.LL4)

  x <- max(drug.col.par[2], drug.row.par[2]) + 1

  for (i in 1:(nrow(drug.col) - 1)) {
    for (j in 1:(nrow(drug.row) - 1)) {
      x1 <- drug.col$dose[i + 1]
      x2 <- drug.row$dose[j + 1]

      options(warn = -1)
      slv <- tryCatch({
        slv <- nleqslv::nleqslv(x, eq, method = "Newton", x1=x1, x2=x2,
                                drug.col.par = drug.col.par,
                                drug.row.par = drug.row.par)
        },
        error = function(){
          slv <- list(termcd = 999)
        }
      )

      if (slv$termcd < 3) {
        y.loewe <- slv$x
      } else {
        y.loewe1 <- fun(x1, x2, drug.col.par, drug.col.fct) # x1 col, x2 row
        y.loewe2 <- fun(x1, x2, drug.row.par, drug.row.fct) # x1 col, x2 row
        y.loewe <- max(y.loewe1, y.loewe2)
      }

      if (y.loewe > 100) {
        y.loewe <- 100
      }

      loewe.mat[j + 1, i + 1] <- y.loewe
    }
  }

  synergy.mat <- response.mat - loewe.mat

  return(synergy.mat)

  options(warn = 0)
  # clean up
  gc()
}