node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
cell <- GetCell(node, "U87", "name")
info <- GetCellInfo(cell, "accession")

ref.list <- XML::xmlChildren(cell[[1]])$`xref-list`
ref <- sapply(XML::xmlChildren(ref.list), XML::xmlAttrs)
ccle <- ref[3, which(ref[1,] == "CCLE")]
tissue <- gsub("_", " ", tolower(gsub("^[^_]+(?=_)_","",ccle, perl = TRUE)))
