# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

#' Generate "cell_line" and "tissue" file
#'
#' @param ids
#' @param type
GenerateCell <- function(name){
  # Get all cell lines' cellosaurus accession
  doc <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
  cell <- GetCell(doc, ids = name, type = "name")
  acc <- GetCellInfo(cell, "accession")

  # Checking whether cell lines have been in DrugComb
  exist.cell <- CheckCell(acc)

  if (length(exist.cell$new) == 0) {
  cell_line <- data.frame(name = character(), synonyms = character(),
                        cellosaurus_accession = character(),
                        disease_id = character(),
                        id = numeric(),
                        tissue_id = numeric())
  disease <- data.frame(name = character(), id = character())
  tissue <- data.frame(id = numeric(), name = character())
  } else {
    cell <- GetCell(doc, ids = as.character(exist.cell$new), type = "accession")
    name <- GetCellInfo(cell, "name")
    dis <- GetCellInfo(cell, "disease")
    dis.uni <- unique(dis)
    tis <- GetCellInfo(cell, "tissue")
    tis.uni <- unique(tis)
    accession <- GetCellInfo(cell, "accession")
    id <- seq(exist.cell$n + 1, length.out = length(exist.cell$new))

    exist.disease <- CheckDisease(as.character(dis.uni[, 2]))
    disease <- exist.disease$new

    exist.tissue <- CheckTissue(tis.uni[,1])
    if (length(exist.tissue$new) == 0) {
      tissue <- data.frame(id = numeric(), name = character())
      tissue_id <- exist.tissue$old[match(tis, exist.tissue$old[ ,1]), 2]
    } else {
      tissue <- data.frame(id = seq(exist.tissue$n + 1,
                           length.out = length(exist.tissue$new)),
                           name = exist.tissue$new)
      tissue_id <- rbind(exist.tissue$old, tissue)
      tissue_id <- tissue_id$id[match(tis, tissue_id$name)]
    }

    cell_line <- cbind(name, accession, dis[, "disease_id"], id, tissue_id)
    colnames(cell_line) <- c("name", "synonyms", "cellosaurus_accession",
                             "disease_id", "id", "tissue_id")
  }
  return(list(cell_line, tissue, disease))
}
