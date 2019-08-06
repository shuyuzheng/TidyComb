# TidyComb
# Functions for calculating cell line's sensitivity to drugs or drug combination
# Copyright Shuyu Zheng
#
# Functions on this page:
#
# Exported:
#
# CalculateSens: Calculate sensitivity score (relative inhibition)
# ImputeIC50: Impute missing value at IC50 concentration of drug
# PredictResponse: Predict response value at certain drug dose
#
# Internal:
# scoreCurve/scoreCurve.L4: facility functions for CalculateSens
# own_log/own_log2: facility functions for CalculateSens
# CalculateIC50: Transform IC50 from coefficients from fitted dose-response model

#' Calculate sensitivity score (relative inhibition)
#'
#' Function \code{CalculateSens} calculates cell line sensitivity to a drug or a
#' combination of drugs from dose response curve.
#'
#' This function measures the sensitivity by calculating the Area Under Curve
#' (AUC) according to certurn dose response curve. The lower bouder is chosen at
#' the first non-zero concentration in the dose response data.
#'
#' @param df A data frame. It contains two variables:
#' \itemize{
#'   \item \strong{dose} the concentrations of drugs.
#'   \item \strong{response} the response of cell lines at crresponding doses.
#'   We use inhibition rate of cell line growth to measure the response.
#' }
#' \strong{Note}: The input data frame must be sorted by "dose" with ascending
#' order.
#'
#' @return A number. The sensitivity score calculated from input dose-response
#' data
#'
#' @export
#'
#' @examples
#' df <- data.frame(dose = c(0, 0.1954, 0.7812, 3.125, 12.5, 50),
#'                  response = c(2.95, 3.76, 18.13, 28.69, 46.66, 58.82))
#' sensitivity <- CalculateSens(df)

CalculateSens <- function(df) {
  #options(show.error.messages = FALSE)
  df <- df[which(df$dose != 0),]
  if (nrow(df) == 1) {
    score <- df$response[1]
  } else {
    score <- tryCatch({
      # Skip zero conc, drc::LL.4()
      # fitcoefs <- drc::drm(formula = as.numeric(df[2:nrow(df),1]) ~
      #                       as.numeric(rownames(df)[2:nrow(df)]),
      #                     fct = drc::LL.4())$coefficients
      # If fit works, call scoreCurve()
      fitcoefs <- drc::drm(response ~ dose, data = df, fct = drc::LL.4(),
                           control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                               otrace = TRUE))$coefficients
      score <- round(scoreCurve(d = fitcoefs[3] / 100,
                                c = fitcoefs[2] / 100,
                                b = fitcoefs[1],
                                m = log10(fitcoefs[4]),
                                c1 = log10(min(df$dose)),
                                c2 = log10(max(df$dose)),
                                t = 0), 3)
    }, error = function(e) {
      # Skip zero conc, log, drc::L.4()
      # message(e)
      fitcoefs <- drc::drm(response ~ log10(dose), data = df, fct = drc::L.4(),
                           control = drc::drmc(errorm = FALSE, noMessage = TRUE,
                                               otrace = FALSE))$coefficients
      score <- round(scoreCurve.L4(d = fitcoefs[3] / 100,
                                   c = fitcoefs[2] / 100,
                                   b = fitcoefs[1],
                                   e = fitcoefs[4],
                                   c1 = log10(min(df$dose)),
                                   c2 = log10(max(df$dose)),
                                   t = 0), 3)
    }
    )
  }
  #options(show.error.messages = TRUE)
  return (score)

  #Clean up
  gc()
}

# CSS - scoreCurve
# New function used to score sensitivities given either a single-agent or a
# fixed conc (combination) columns.
# The function calculates the AUC of the log10-scaled dose-response curve.

# IMPORTANT: note that with L.4() calls, this value is already logged since the
# input concentrations are also logged.
# c1:log10(min conc) (this is the minimal nonzero concentration)
# c2:log10(max conc) (this is the maximal concentration)
# t: threshold (usually set to zero)
# own_log: log(1+10^(b*(c-x)))

scoreCurve <- function(b, c, d, m, c1, c2, t) {
  # y <- c + (d - c) / (1 + (e / x) ^ (-b)) # LL.4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)
  int_y <- (((((d - c) * own_log(-b, c2, m)) / ((-b) * log(10))) + c * c2) -
              ((((d - c) * own_log(-b, c1, m)) / ((-b) * log(10))) + c * c1))
  # int_y <- (((((a - d) * own_log(b, c, x2)) / (b * log(10))) + a * x2) -
  #            ((((a - d) * own_log(b, c, x1)) / (b * log(10))) + a * x1)) -
  #              (t * (x2 - x1))

  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}

own_log = function(b, c, x)
{
  arg = 1 + 10^(b*(c-x))
  if(is.infinite(arg)==T) res = b*(c-x)*log(10) else res = log(arg)
  return(res)
}
# CSS - scoreCurve for L.4() model
# New function used to score sensitivities given either a single-agent or a
# fixed conc (combination) columns.
# The function calculates the AUC of the log10-scaled dose-response curve
# own_log2: log(1+exp(x))

scoreCurve.L4 <- function(b, c, d, e, c1, c2, t) {
  # y <- c + (d - c) / (1 + exp(b * (x - e))) # L4
  # b <- coef[1]
  # c <- coef[2]
  # d <- coef[3]
  # e <- coef[4]
  # m <- log10(e)

  int_y <- d * (c2 - c1) + ((c - d) / b) *
    (own_log2(b * (c2 - e)) - own_log2(b * (c1 - e)))
  ratio <- int_y / ((1 - t) * (c2 - c1))
  sens <- ratio * 100 # scale by 100
  return(sens)
}

# own_log2, natural base
# calculate log(1+exp(x)) to be used in scoreCurve.L4 function
own_log2 = function(x)
{
  arg = 1 + exp(x)
  if(is.infinite(arg)==T) res = x else res = log(arg)
  return(res)
}

#' Impute missing value at IC50 concentration of drug
#'
#' \code{ImputeIC50} uses the particular experiment's values to predict the
#' missing values at the desired IC50 concentration of the drug.
#
#' This function is only called when trying to fix a drug at its selected IC50
#' concentration where the response values have not been tested in experiment.
#'
#' \code{ImputeIC50} fits dose-response models (with \code{\link[drc]{drm}}
#' function) by fixing the concentrations of the
#' \strong{other} drug successively, and uses each fit to predict the missing
#' value at the combination (missing IC50, fixed conc).
#'
#' @param response.mat A matrix. It contains response value of a block of drug
#' combination.
#'
#' @param row.ic50 A numeric. The IC50 value of drug added to rows.
#'
#' @param col.ic50 A numeric. The IC50 value of drug added to columns.
#'
#' @return a data frame contains all response value at the IC50 concentration
#' of certein drug. It could be directly passed to function
#' \code{CalculateSens} for scoring.
#'
#' @export
ImputeIC50 <- function(response.mat, col.ic50, row.ic50) {

  colconc <- as.numeric(colnames(response.mat))
  rowconc <- as.numeric(rownames(response.mat))
  n_col <- length(colconc)
  n_row <- length(rowconc)

  if (n_row == 2) {
    tempcf_c <- data.frame(dose = colconc, response = response.mat[2, ])
  } else {
    response <- apply(response.mat, 2, function(x){
          df <- data.frame(dose = rowconc, response = x)
          pred <- PredictResponse(df, row.ic50)
          return(pred)
        }
      )
    tempcf_c <- data.frame(dose = colconc, response = response)
  }

  if (n_col == 2) {
    tempcf_r <- data.frame(dose = rowconc, response = response.mat[, 2])
  } else {
    response <- apply(response.mat, 1, function(x){
        df <- data.frame(dose = colconc, response = x)
        pred <- PredictResponse(df, col.ic50)
        return(pred)
      }
    )
    tempcf_r <- data.frame(dose = rowconc, response = response)
  }

  tempres <- list(tempcf_c = tempcf_c, tempcf_r = tempcf_r)
  return(tempres)

  # Clean up
  gc()
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

#' Predict response value at certain drug dose
#'
#' \code{PredictResponse} uses \code{\link[drc]{drm}} function to fit the dose
#' response model and generate the predict response value at the given dose.
#'
#' \strong{Note}: Random number generator used in \code{AddNoise} with
#' \code{method = "random"}. If the analysis requires for reproductiblity,
#' plesase set the random seed before calling this function.
#'
#' @param df A data frame. It contains two variable:
#' \itemize{
#'   \item \strong{dose} a serial of concentration of drug;
#'   \item \strong{response} the cell line response to each concentration of
#'   drug. It should be the inhibition rate according to negative control.
#' }
#'
#' @param dose A numeric value. It specifies the dose at which user want to
#' predict the response of cell line to the drug.
#'
#' @return A numeric value. It is the response value of cell line to the drug at
#' inputted dose.
#'
#' @author Shuyu Zheng{shuyu.zheng@helsinki.fi}
#'
#' @export
PredictResponse <- function(df, dose) {
  if (stats::var(df$response, na.rm = TRUE) == 0) {
    pred <- df$response[1]
  } else {
    model <- synergyfinder::FitDoseResponse(df)

    if (model$call$fct[[1]][[3]] == "LL.4") {
      pred <- stats::predict(model, data.frame(dose = dose))
    } else {
      pred <- stats::predict(model, data.frame(dose = log(dose)))# NB! use log
    }

    if (pred > 100) {
      pred <- 100 + stats::runif(1, -0.01, 0)
    }
  }
  return(pred)
}