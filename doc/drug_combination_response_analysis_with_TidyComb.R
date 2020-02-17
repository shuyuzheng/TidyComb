## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TidyComb)
library(reshape2)
library(dplyr)

## ----message=FALSE, warning=FALSE---------------------------------------------
df <- data.frame(dose = c(0, 0.1954, 0.7812, 3.125, 12.5, 50),
                response = c(2.95, 3.76, 18.13, 28.69, 46.66, 58.82))
sens <- CalculateSens(df)
print(sens)

## ----message=FALSE, warning=FALSE---------------------------------------------
pred <- CalculateSens(df, pred = TRUE)$pred
RIConfidenceInterval(pred)

## -----------------------------------------------------------------------------
data <- read.csv(system.file("template.csv", package = "TidyComb"),
                 stringsAsFactors = FALSE)
response.mat <- reshape2::acast(conc_r~conc_c, value.var = "inhibition",
                               data = data[data$block_id == 1, ])
response.mat

## ----message=FALSE, warning=FALSE---------------------------------------------
CalculateMat(response.mat, summary.only = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
res <- CalculateMat(response.mat)
res

## -----------------------------------------------------------------------------
data

## ----message=FALSE------------------------------------------------------------
CalculateTemplate(data, summary.only = TRUE)

## -----------------------------------------------------------------------------
res <- ParCalculateTemplate(data, cores = 4)
res

