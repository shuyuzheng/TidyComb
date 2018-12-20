# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng


# Get all cell lines' cellosaurus accession
MatchCellAcc <- function(name){
  doc <- GetAllCell()
  cell <- GetCell(doc, ids = name, type = "name")
  acc <- GetCellInfo(cell, "accession")
  names <- GetCellInfo(cell, "name")
  all_name <- NA
  for (i in 1:nrow(names)) {
    all_name[i] <- paste(names[i, 1], names[i, 2], sep = "; ")
  }

  info <- data.frame(acc, all_name, stringsAsFactors = FALSE)
  df <- NULL
  for (i in name) {
    temp <- NULL
    temp <- info[grep(i, info[, "all_name"]), ]
    if (nrow(temp) == 0) {
      temp[1, ] = c(NA, NA)
    }
    temp$input_name <- i
    df <- rbind(df, temp)
  }
  df <- df[, c("input_name", "cellosaurus_accession", "all_name")]
  message("Please check the table and make sure that each cell line matches",
          " with correct cellosaurus accessions. \n",
          "Query on https://web.expasy.org/cellosaurus/ may be helpful.")
  return(df)
}



#' Generate "cell_line" and "tissue" file
#'
#' @param acc A vector of charactors contains cellosaurus accessions of
#'  interested cell lines.
#'
#' @return  A list with 3 data frame contains
#' \itemize{
#'   \item \strong{cell_line} The cell_line table prepared for uploading.
#'   \item \strong{tissue} The tissue table prepared for uploading.
#'   \item \strong{disease} The disease table prepared for uploading.
#' }
#' @export
#'
GenerateCell <- function(acc){
  message("Generating cell_line, disease, tissue...")
  if (!all(grepl("CVCL_.+", acc))) {
    stop("Illegal cellosaurus accession! Please check your input list. \n\r",
         "Cellosaurus accessions should in form 'CVCL_xxxx'.")
  }
  doc <- GetAllCell()
  # Checking whether cell lines have been in DrugComb
  exist.cell <- CheckCell(acc)

  if (length(exist.cell$new) == 0) {
  cell_line <- data.frame(name = character(),
                          synonyms = character(),
                          cellosaurus_accession = character(),
                          disease_id = character(),
                          id = numeric(),
                          tissue_id = numeric(),
                          stringsAsFactors = FALSE)
  disease <- data.frame(name = character(),
                        id = character(),
                        stringsAsFactors = FALSE)
  tissue <- data.frame(id = numeric(),
                       name = character(),
                       stringsAsFactors = FALSE)
  } else {
    cell <- GetCell(doc, ids = as.character(exist.cell$new), type = "accession")
    names <- GetCellInfo(cell, "name")
    dis <- GetCellInfo(cell, "disease")
    dis.uni <- unique(dis)
    tis <- GetCellInfo(cell, "tissue")
    tis.uni <- stats::na.omit(unique(tis))
    accession <- GetCellInfo(cell, "accession")
    id <- seq(exist.cell$n + 1, length.out = length(exist.cell$new))

    exist.disease <- CheckDisease(as.character(dis.uni[, 2]))
    disease <- exist.disease$new

    exist.tissue <- CheckTissue(tis.uni[,1])
    if (length(exist.tissue$new) == 0) {
      tissue <- data.frame(id = numeric(),
                           name = character(),
                           stringsAsFactors = FALSE)
      tissue_id <- exist.tissue$old[match(tis, exist.tissue$old[ ,1]), 2]
    } else {
      tissue <- data.frame(id = seq(exist.tissue$n + 1,
                                    length.out = length(exist.tissue$new)),
                           name = exist.tissue$new,
                           stringsAsFactors = FALSE)
      tissue_id <- rbind.data.frame(exist.tissue$old, tissue)
      tissue_id <- tissue_id$id[match(tis, tissue_id$name)]
    }

    cell_line <- data.frame(names, accession, dis[, "disease_id"],
                                  id, tissue_id, stringsAsFactors = FALSE)
    colnames(cell_line) <- c("name", "synonyms", "cellosaurus_accession",
                             "disease_id", "id", "tissue_id")
  }

  # Generate cell_id index table
  cell_id <- cell_line[, c("cellosaurus_accession", "id")]
  cell_id <- rbind.data.frame(cell_id, exist.cell$old)
  return(list(cell_line = cell_line,
              tissue = tissue,
              disease = disease,
              cell_id = cell_id))
}
