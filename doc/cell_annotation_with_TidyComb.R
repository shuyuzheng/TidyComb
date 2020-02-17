## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TidyComb)

## -----------------------------------------------------------------------------
CellVersion()

## -----------------------------------------------------------------------------
CellVersion(file = system.file("cellosaurus.xml", package = "TidyComb"))

## -----------------------------------------------------------------------------
cell_name <- c("U87", "HSTS", "LNCAP")
MatchCellAcc(cell_name, file = system.file("cellosaurus.xml", package = "TidyComb")) # I'm using the embedded data. You can change is into the file downloaded by yourself

## -----------------------------------------------------------------------------
doc <- ParseCell(file = system.file("cellosaurus.xml", package = "TidyComb")) # I'm using the embedded data. You can change is into the file downloaded by yourself
GetCellInfo(c("CVCL_0022", "CVCL_0395"), doc)

