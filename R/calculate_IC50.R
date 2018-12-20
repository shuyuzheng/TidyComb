# CSS - calculateIC50
# This function takes in a data frame as input
# The data frame must contain one column of response (inhibition values)
# The row names of the data frame are the increasing concentrations of the
# tested drug The function tries fitting a drm model using one of three
# consecutive methods dertailed below as soon as one of the models converges,
# the IC50 parameter is returned

calculateIC50 <- function(df) {
  # options(show.error.messages = F)
  # df$conc[1] = 10^-10
  ic50 <- tryCatch({
      # drm(response ~ conc, data = df, fct = drc::LL.4())$coefficients[4]
      res <- drc::drm(response ~ conc, data = df, fct = drc::LL.4())
      # b <- res$coefficients[1]
      # c <- res$coefficients[2]
      # d <- res$coefficients[3]
      e <- res$coefficients[4]
      ic50 <- e
    }, error = function(e){
      df$conc[1] <- 10^-10
      # ic50 <- drc::drm(response ~ log(conc), data = df,
      #                  fct = drc::L.4())$coefficients[4]
      res <- drc::drm(response ~ log(conc), data = df, fct = drc::L.4())
      # b <- res$coefficients[1]
      # c <- res$coefficients[2]
      # d <- res$coefficients[3]
      e <- res$coefficients[4]
      # remember to take it back to natural scale
      ic50 <- exp(e)
    }
  )
  # options(show.error.messages = T)
  if (ic50 > max(df$conc)) ic50 = max(df$conc)
  return (ic50)
}
