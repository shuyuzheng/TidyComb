# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng


#' Check the Cellosaurus data version.
#'
#' Cellosaurus publishes its dataset as an XML file which could be downloaded
#' from ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml. This
#' could be used to check the publish version of this file (by default) or a
#' XML downloaded before (by passing the "file" argument).
#'
#' @param file A character. The xml file contains the Cellosaurus dataset.
#' @examples
#' # Check the online Cellosaurus database version.
#' CellVersion()
#'
#' # Check the local Cellosaurus XML document version.
#' CellVersion(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
CellVersion <- function(
          file = "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml",
          ...){
  #build handler
  startElement = function(ctxt, name, attrs, ...) {
    if (name == "release") {
      print(attrs)
      XML::xmlStopParser(ctxt)
    }
  }
  class(startElement) = "XMLParserContextFunction"

  print("Cellosaurus database version is:")
  invisible(XML::xmlEventParse(file,
                               handlers = list(startElement = startElement)))
}


#' Loading a Cellosaurus XML dataset
#'
#' Cellosaus published "cellosaurus.xml" contains 5 part of information:
#' "header", "cell-line-list", "publication-list", "copyright". (more
#' information in "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xsd")
#' This function will parse the XML file and extract "cell-line-list" node which
#' includes all information about cell lines as a XML document object.
#'
#' @param file A Cellosaurus XML file. It could be the URL point to online file
#' ("ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml") or a local
#' file already downloaded (recommend).
#'
#' @return An XMLNode containing all cell lines' information archieved in
#' Cellosaurus dataset.
#'
#' @example
#' CellDoc <- GetAllCell(system.file("extdata",
#'                                   "cellosaurus.xml",
#'                                   package = "TidyComb"))
#'
#' @export
GetAllCell <- function(file){
    doc <- XML::xmlInternalTreeParse(file)
    all.cells <- XML::xmlRoot(doc)[[2]]
    return(all.cells)
}


#' Find matitching cell-lines
#'
#' This function searches the value of <name> tags in Cellosaurus XML document
#' to find entries that matches the names provided in the \code{names}
#' parameter.
#'
#' @param node An XMLNodelist. It is the output of \code{\link{GetAllCell}}
#' function which contains all cell lines' information in Cellosaurus dataset.
#'
#' @param names A character or a vector of characters. It is the name of cell
#' lines that will be searched in Cellosaurus XML file.
#'
#' @return An XMLNodeSet containing matched cell lines in the dataset. If no
#' cell line is matched, a NULL list will be return.
#'
#' @examples
#' all.cells <- GetAllCell(system.file("extdata",
#'                                     "cellosaurus.xml",
#'                                      package = "TidyComb"))
#' cell.lines <- GetCell(all.cells, c("U87", "A549"))
#'
#' @export
GetCell <- function(node, names){
      xpath <- paste0("//name[text()='",
                    paste(names, collapse = "' or text() = '"),
                    "']/ancestor::cell-line")
    cells <- XML::getNodeSet(node, xpath)
}

#' Extract cell line information
#'
#' This function will extract the information required by DrugComb database (
#' including cell name, cell synonyms, disease, disease_id ) of the cell lines.
#' The input should be "XMLInternalElementNode" that contains the information
#' of the cell lines. It could be output of either  \code{\link{GetAllCell}} (
#' all Cellosaurus archived cell lines) or \code{\link{GetCell}} (cell lines
#' you're interested in).
#'
#' @param node An "XMLInternalElementNode" extracted from the Cellosaurus XML
#' file by either \code{\link{GetAllCell}} or \code{\link{GetCell}}
#'
#' @param term A vector of character. The type of information you want to
#' extract. A valible values are: "synonyms", "accession", "disease",
#' "diseaseID", "name", "all".
#'
#' @return A data frame containing the information selected by \code{term}
#' extracted from the \code{node} item.
#'
#' @examples
#' all.cells <- GetAllCell(system.file("extdata",
#'                                     "cellosaurus.xml",
#'                                      package = "TidyComb"))
#' cell.lines <- GetCell(all.cells, c("U87", "A549"))
#' cell.info <- GetCellInfo(cell.lines, "synonyms")

GetCellInfo <- function(node, term) {

  terms <- c("synonyms", "accession", "disease", "diseaseID", "name", "all")

  if (!(term %in% terms)) {
    stop("Can not find term '", term, "'. Please chose one term from: ",
         paste(terms, collapse = ", "), ".")
  }

  if (term == "all") {
    xpath <- paste(paste0("//", terms[-length(terms)]), collapse = " | ")
  }





  info <- NULL
  stepi <- 0
  n <- length(cells)

  for (cell in cells) {
    # set a pregressing index
    message(round(stepi/n * 100), "%", " is completed",
            "\r", appendLF = FALSE)
    flush.console()

    # building containers
    synonyms <- character()
    disease <- character()
    disease_id <- character()
    accession <- character()
    name <- character()

    cell <- XML::xpathApply(node, path = paste0("//name[text()=\'", names,
                                          "\']/../.."),xmlToList)

    node.info <- NULL
    for (i in length(node)) {
      new_item <- NA
      # ger synonyms
      new_item <- unlist(sapply(node[[i]]$`name-list`, "[", "text"))
        if (is.null(new_item)) {
            new_item <- NA
          }else{
            new_item <- paste(new_item[-1], collapse = ";")
          }

      synonyms <- c(synonyms, new_item)
      new_item <- NA

      # get accession ID
      new_item <- unlist(sapply(node[[i]]$"accession-list","[", "text"))
        if (is.null(new_item)) {
         new_item <- NA
        }
      accession <- c(accession, new_item)
      new_item <- NA

      # get disease type
      new_item <- unlist(node[[i]]$"disease-list")
      if (length(new_item) == 0) {
        new_item <- NA
      }
      disease <- c(disease, new_item["cv-term.text"])
      disease_id <- c(disease_id, new_item["cv-term..attrs.accession"])
      new_item <- NA

      # record the cell name
      name <- c(name, cell)

      # combine the data of one node
      node.info <- data.frame(name, synonyms, cellosaurus_accession = accession,
                         disease, disease_id, stringsAsFactors = FALSE)
    }

    # combine the data of one cell
    info <- rbind.data.frame(info, node.info)

    stepi <- stepi + 1
    }

  # clean

  gc()
  return(info)
}
