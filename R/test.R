node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
cell <- GetCell(node, "U87", "name")
info <- GetCellInfo(cell, "tissue")

ref.list <- XML::xmlChildren(cell[[1]])$`xref-list`
ref <- sapply(XML::xmlChildren(ref.list), XML::xmlAttrs)
ccle <- ref[3, which(ref[1,] == "CCLE")]
tissue <- gsub("_", " ", tolower(gsub("^[^_]+(?=_)_","",ccle, perl = TRUE)))

drug <- read.csv(system.file("extdata", "drug.csv", package = "TidyComb"),
                 header = FALSE,
                 col.names = c("id", "cid"))

template <- read.csv(system.file("extdata", "template.csv", package = "TidyComb"),
                     stringsAsFactors = FALSE)
cids <- unique(c(template$drug_row_cid, template$drug_col_cid)) %>% na.omit()
