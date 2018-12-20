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
# for scoring

imputeData2 <- function(response.mat, ic_c, ic_r) {

  colconc <- as.numeric(colnames(response.mat))
  rowconc <- as.numeric(rownames(response.mat))
  m_col <- length(colconc)
  m_row <- length(rowconc)

  if (m_row == 2) {
    tempcf_c <- data.frame(conc = as.numeric(colnames(response.mat)),
                           response = response.mat[m_row, ])
  } else {
  response <- apply(response.mat, 2, function(x) {
    data <- data.frame(conc = as.numeric(rownames(response.mat)), response = x)
    if (stats::var(x, na.rm = TRUE) == 0) {
      pred <- x[1]
    } else {
      pred <- tryCatch({
          # model <- drc::drm(response ~ conc,
          #                   data = data[-which(data$conc == 0), ],
          #                   fct = drc::LL.4())
          model <- drc::drm(response ~ conc, data = data, fct = drc::LL.4())

          stats::predict(model, data.frame(conc = ic_r))
        }, error = function(e){
          # model <- drc::drm(response ~ log(conc),
          #              data = data[-which(data$conc == 0), ],
          #              fct = drc::L.4())
          data$conc[1] = 10^-10
          model = drc::drm(response ~ log(conc), data = data, fct = drc::L.4())

          stats::predict(model, data.frame(conc = log(ic_r))) # NB! use log
        }
      )
      if (pred > 100) {
        pred <- 100 + stats::runif(1, -0.01, 0)
      }
    }
    return(pred)
   }
  )
  tempcf_c <- data.frame(conc = as.numeric(colnames(response.mat)),
                         response = response)
  }

  if (m_col == 2) {
    tempcf_r <- data.frame(conc = as.numeric(rownames(response.mat)),
                           response = response.mat[, m_col])
  } else {
    response <- apply(response.mat, 1, function(x) {
      data <- data.frame(conc = as.numeric(colnames(response.mat)),
                         response = x)
      if(stats::var(x, na.rm = T) == 0) {
        pred = x[1]
      } else {
        pred <- tryCatch({
            # model <- drc::drm(response ~ conc,
            #                   data = data[-which(data$conc == 0), ],
            #                   fct = drc::LL.4())
            model <- drc::drm(response ~ conc, data = data, fct = drc::LL.4())
            stats::predict(model, data.frame(conc = ic_c))
          }, error = function(e) {
            data$conc[1] <- 10^-10
            model <- drc::drm(response ~ log(conc), data = data,
                              fct = drc::L.4())
            # model <- drc::drm(response ~ log(conc),
            #                   data = data[-which(data$conc == 0), ],
            #                   fct = drc::L.4())
            stats::predict(model, data.frame(conc = log(ic_c))) # NB! use log
          }
        )
        if (pred > 100) {
          pred <- 100 + stats::runif(1, -0.01, 0)
        }
      }
      return(pred)
    }
    )
    tempcf_r <- data.frame(conc = as.numeric(rownames(response.mat)),
                           response = response)
  }
  tempres <- list(tempcf_c = tempcf_c, tempcf_r = tempcf_r)
  return(tempres)
}
