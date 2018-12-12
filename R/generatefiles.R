# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

#' Generate "cell_line" and "tissue" file
#'
#' @param ids
#' @param type
GenerateCell <- function(name){
  # Get all cell lines' cellosaurus accession
  doc <- GetAllCell(system.file("extdate", "database", "cellosaurus.xml", package = "TidyComb"))
  cell <- GetCell(doc, ids = name, type = "name")
  acc <- GetCellInfo(cell, "accession")

  # Checking whether cell lines have been in DrugComb
  exist <- CheckCell(accession = acc)

  if (length(exist$new) == 0) {
  cell_line <- data.frame(name = character(), synonyms = character(),
                        cellosaurus_accession = character(),
                        disease_id = character(),
                        id = numeric(),
                        tissue_id = numeric())
  disease <- data.frame(name = character(), id = character())
  tissue <- data.frame(id = numeric(), name = character())
  } else {
    cell <- GetCell(doc, ids = exist$new, type = "accession")
    name <- GetCellInfo(cell, "name")
    dis <- GetCellInfo(cell, "disease")
    tis <- GetCellInfo(cell, "tissue")
    accession <- GetCellInfo(cell, "accession")
    id <- seq(exist$n + 1, length.out = length(exist$new))

    exist.disease <- CheckDisease(dis)
    disease <- exist.disease$new

    exist.tissue <- CheckTissue(tis)
    tissue <- data.frame(id = seq(exist.tissue$n + 1,
                                  length.out = length(exist.tissue$new)),
                         name = exist.tissue$new)
    tissue_id <- cbind.data.frame(exist.tissue$old, tissue$id)
    tissue_id <- tissue_id$id[match(tis, tissue_id$name)]

    cell_line <- cbind.data.frame(name, accession, dis$disease_id,
                                  id, tissue_id)
    colnames(cell_line) <- c("name", "synonyms", "cellosaurus_accession",
                             "disease_id", "id", "tissue_id")
  }
  return(list(cell_line, tissue, disease))
}
