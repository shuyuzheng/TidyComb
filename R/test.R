node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
cell <- GetCell(node, "name", c("U251"))
info <- GetCellInfo(cell)

accession <- XML::xmlValue(XML::xmlChildren(cell[[1]])$`accession-list`)

