# CSS - imputeData
# This function is only called when trying to fix a drug at its selected IC50
# concentration in an experiment where it hasn't been tested with that
# concentration.
# The purpose of the function is then to use the particular experiment's values
# to predict the missing values at the desired IC50 concentration of the drug.
# The function fits drc models by fixing the concentrations of the OTHER drug
# successively, and uses each fit to predict the missing value at the
# combination (missing IC50, fixed conc).
# It finally returns a data frame ready to be passed to computeSensitivity()
# for scoring.

imputeData <- function(tempmaindf, dimension_to_fix,
                       dimension_to_predict, dose_to_predict) {

  # get the concentration range
  concset <- unique(tempmaindf[,
                               which(colnames(tempmaindf) == dimension_to_fix)])

  # create vector for predicted values
  pred_vals <- c()

  # fit drm functions
  m <- length(concset)
  for (i in 1:m) {
    tempdf <- tempmaindf[which(tempmaindf[,
                               which(colnames(tempmaindf) == dimension_to_fix)]
                                     == concset[i]),
                         c(dimension_to_predict,'response')]
    names(tempdf)[1] <- "conc"
    # options(show.error.messages = FALSE)
    if (nrow(tempdf) == 2) {
      newval = tempdf$response[2]
    } else {
      newval <- tryCatch({
          # Skip zero conc, drc::LL.4()
          # modelx-drc::drm(formula = as.numeric(tempdf[2:nrow(tempdf),1]) ~
          #                 as.numeric(rownames(tempdf)[2:nrow(tempdf)]),
          #                 fct = drc::LL.4())
          model <- drc::drm(response ~ conc,
                            data = tempdf[-which(tempdf$conc == 0), ],
                            fct = drc::LL.4()) # remove the zero concentration
          newval <- stats::predict(model, data.frame(conc = dose_to_predict))
        }, error <- function(e) {
          # Skip zero conc, log, drc::L.4()
          # model <- drc::drm(formula = as.numeric(tempdf[2:nrow(tempdf),1]) ~
          #                   log(as.numeric(rownames(tempdf)[2:nrow(tempdf)])),
          #                   fct = drc::L.4())
          model <- drc::drm(response ~ log(conc),
                            data = tempdf[-which(tempdf$conc == 0), ],
                            fct = drc::L.4()) # remove the zero concentration
          newval <- stats::predict(model,
                                   data.frame(conc = log(dose_to_predict)))
        }
      )
    }
    # options(show.error.messages = T)
    if (newval > 100) {
      newval <- 100 + stats::runif(1, -0.01, 0.01)
    } # random noise to avoid singularity
    pred_vals <- c(pred_vals, newval)
  }

  tempres <- data.frame(conc = concset, response = pred_vals)
  return(tempres)
}
