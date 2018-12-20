# CSS - calculateIC50
# This function takes in a data frame as input
# The data frame must contain one column of response (inhibition values)
# The row names of the data frame are the increasing concentrations of the tested drug
# The function tries fitting a drm model using one of three consecutive methods dertailed below
# As soon as one of the models converges, the IC50 parameter is returned

CalculateAbsoluteIC50 <- function(df) {
  # options(show.error.messages = F)
  df$conc[1] <- 10^-10
  ic50 <- tryCatch({
      #No cutoff, drc::LL.4()
      # as.numeric(drc::drm(formula = as.numeric(df$response) ~
      #                     as.numeric(rownames(df)),
      #                     fct = drc::LL.4())$coefficients[4])
      model <- drc::drm(response ~ conc, data = df, fct = drc::LL.4())
      res <- LL4Reverse(model$coefficients[1], model$coefficients[2],
                        model$coefficients[3], model$coefficients[4], 50)
      # res <- ED(model, 50, type = "absolute")[1]
      max.conc <- max(df$conc)
      # if a drug cannot reach 50% then take the maximal tested concentration
      ifelse(is.nan(res) | res > max.conc, max.conc, res)
    }, error <- function(e) {
      # ic50 <- tryCatch({
      #     cutoff, drc::LL.4()
      #     as.numeric(drc::drm(formula = as.numeric(df$response) ~
      #                         c(10^-10, as.numeric(rownames(df)[2:nrow(df)])),
      #                         fct = drc::LL.4())$coefficients[4])
      #    df$conc[1] <- 10^-10
      #    drc::drm(response~conc, data = df, fct = drc::LL.4())$coefficients[4]
      #    }, error = function(e1) {
      #      cutoff, log, drc::L.4()
      #     ic50 <- as.numeric(drc::drm(formula = as.numeric(df$response) ~
      #                               c(log(10^-10),
      #                             log(as.numeric(rownames(df)[2:nrow(df)]))),
      #                             fct = drc::L.4())$coefficients[4])
      #     df$conc[1] <- 10^-10
      model <- drc::drm(response ~ log(conc), data = df, fct = drc::L.4())
      res <- L4Reverse(model$coefficients[1], model$coefficients[2],
                      model$coefficients[3], model$coefficients[4], 50)
      # ED function has a bug, no difference between absolute and relative types
      # res <- ED(model, 50, type = "absolute")[1]
      max.conc <- max(df$conc)
      exp.res <- exp(res)
      # if a drug cannot reach 50% then take the maximal tested concentration
      ifelse(is.nan(res) | exp.res > max.conc, max.conc, exp.res)
      #  }
      # )
    }
  )
  # options(show.error.messages = TRUE)
  return (ic50)
}
