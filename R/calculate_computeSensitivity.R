# CSS - computeSensitivity
# This function takes in a data frame as input
# The data frame has one column containing the response (inhibition) percentages
# The rownames of the data frame are the increasing concentrations of the drug used
# The selected column could be either that of the single-agent
# (in which case this is a DSS calculation)
# Otherwise, it is that in the presence of the second drug fixed at its IC50 score defined earlier
# (in which case this is a CSS calculation)

computeSensitivity <- function(df) {
  #options(show.error.messages = FALSE)
  if (nrow(df) == 2) {
    score <- df$response[2]
  } else {
    score <- tryCatch({
        # Skip zero conc, drc::LL.4()
        # fitcoefs <- drc::drm(formula = as.numeric(df[2:nrow(df),1]) ~
        #                       as.numeric(rownames(df)[2:nrow(df)]),
        #                     fct = drc::LL.4())$coefficients
        # If fit works, call scoreCurve()
        fitcoefs <- drc::drm(response ~ conc, data = df[-which(df$conc == 0),],
                             fct = drc::LL.4())$coefficients
        score <- round(scoreCurve(d = fitcoefs[3] / 100,
                                  c = fitcoefs[2] / 100,
                                  b = fitcoefs[1],
                                  m = log10(fitcoefs[4]),
                                 c1 = log10(df$conc[2]),
                                 c2 = log10(df$conc[nrow(df)]),
                                  t = 0), 3)
      }, error = function(e) {
        # Skip zero conc, log, drc::L.4()
        # message(e)
        fitcoefs <- drc::drm(response ~ log10(conc),
                             data = df[-which(df$conc == 0), ],
                             fct = drc::L.4())$coefficients
        score <- round(scoreCurve.L4(d = fitcoefs[3] / 100,
                                     c = fitcoefs[2] / 100,
                                     b = fitcoefs[1],
                                     e = fitcoefs[4],
                                    c1 = log10(df$conc[2]),
                                    c2 = log10(df$conc[nrow(df)]),
                                     t = 0), 3)
      }
    )
  }
  #options(show.error.messages = TRUE)
  return (score)
}
