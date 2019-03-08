# TidyComb
# Functions for fitting dose-response model trough whole response matrix
# Copyright Shuyu Zheng
#
# Functions on this page:
# FitWholeMat

PredictResponse <- function(df, dose) {
  if (stats::var(df$response, na.rm = TRUE) == 0) {
    pred <- df$response[1]
  } else {
    model <- FitDoseResponse(data)

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

FitWholeMat <- function(response.mat) {

}