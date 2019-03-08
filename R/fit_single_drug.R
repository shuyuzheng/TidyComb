# TidyComb
# Functions for fitting single drug dose-response curve
# Copyright Shuyu Zheng
#
# Functions on this page:
#
# FitDoseResponse: Fitting single drug dose-response model

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
  if (nrow(data) != 1 & stats::var(data$response) == 0) {
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

FindModelType <- function(model) {
  type <- model$call$fct[[1]][[3]]
  return(type)
}

CalculateIC50 <- function(coef, type, max.conc){
  if (type == "LL.4") {
    ic50 <- coef[["e:(Intercept)"]]
  } else if (type == "L.4") {
    ic50 <- exp(coef[["e:(Intercept)"]])
  }

  if (ic50 > max.conc) {
    ic50 = max.conc
  }

  return (ic50)

}