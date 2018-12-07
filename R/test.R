main <- function(response){
  cellline <- response$
node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
cell <- GetCell(node, "name", "U87")
info <- GetCellInfo(cell)
return(info)
}

accession <- XML::xmlValue(XML::xmlChildren(cell[[1]])$`accession-list`)

ref.list <- XML::xmlChildren(cell[[1]])$`xref-list`
ref <- sapply(XML::xmlChildren(ref.list), XML::xmlAttrs)
ccle <- ref[3, which(ref[1,] == "CCLE")]
tissue <- str_sub(ccle, grep("_", ccle, fixed = TRUE), str_length(ccle))
gsub("^[^_]+(?=_)","",ccle)

