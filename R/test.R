
cell <- GetCell(GetAllCell(), "U87", "name")
info <- GetCellInfo(cell, "tissue")

id <- "RDHQFKQIGNGIED-UHFFFAOYSA-N"
# calculate_imputeData2 line 41, 70, 77 predict() function
#
input <- read.csv(system.file("extdata",
                              "template.csv",
                               package = "TidyComb"),
                  stringsAsFactors = FALSE)
