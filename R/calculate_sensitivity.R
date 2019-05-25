# TidyComb
# Functions for calculating cell line's sensitivity to drugs or drug combination
# Copyright Shuyu Zheng
#
# Functions on this page:
#

#' Calculate sensitivity score
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
#'
#' @return A number. The sensitivity score calculated from input dose-response
#' data
#'
#' @export

CalculateSens <- function(df) {
  #options(show.error.messages = FALSE)
  df <- df[-which(df$dose == 0),]
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
                                c1 = log10(df$dose[1]),
                                c2 = log10(df$dose[nrow(df)]),
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
                                   c1 = log10(df$dose[1]),
                                   c2 = log10(df$dose[nrow(df)]),
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
