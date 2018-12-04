# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

#' Extract primary name and synonyms of one cell line.
#'
#' This function extracts the primary name (as "name") and synonyms(as"synonyms"
#' ) of one cell line from an XMLIntervalElementNode object extracted from
#' Cellosaurus xml file. Combining this function with \conde{\link[base]{apply}}
#' or \code{\link[base]{sapply}} can extract names information from an
#' XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A data frame contains two variables:
#'   \item{name}{the primary name of cell line.}
#'   \item{synonyms}{synonyms of the cell line separated by semicolons.}
#'
#' @examples
#' node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
#' cell <- GetCell(node, "name", c("U251", "U87"))
#'
#' # get first cell line's name and synonyms
#' name <- GetNames(cell[[1]])
#'
#' # get all cell lines' name and synonyms
#' names <- sapply(cell, GetNames)
GetNames <- function(node) {
  name.list <- XML::xmlToDataFrame(XML::xmlChildren(node)$'name-list',
                              stringsAsFactors = FALSE)
  name <- name.list[1,]
  synonyms <- sapply(name.list,
                     function(x) {
                       paste(x[-1], collapse = "; ")
                     })
  names <- data.frame(name = name, synonyms = synonyms)
  return(names)
}

#' Extract the disease information of one cell line.
#'
#' This function extract the cell line associated disease and disease ID (NCI
#' Thesaurus entry code) which the cell line associated with from an
#' XMLIntervalElementNode object extracted from Cellosaurus xml file. Combining
#' this function with \conde{\link[base]{apply}} or \code{\link[base]{sapply}}
#' can extract disease information from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A data frame contains two variables:
#' \item{disease}{name of the cell line associated disease.}
#' \item{disease_id}{NCI Thesaurus entry code of the associated disease.}
#'
#' @examples
#' node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
#' cell <- GetCell(node, "name", c("U251", "U87"))
#'
#' # get first cell line associated disease and disease ID
#' disease <- GetDisease(cell[[1]])
#'
#' # get all cell lines associated diseases and disease IDs
#' diseases <- sapply(cell, GetDisease)
GetDisease <- function(node){
  disease.list <- XML::xmlChildren(node)$'disease-list'
  disease <- XML::xmlValue(disease.list)
  disease.id <- sapply(XML::xmlChildren(disease.list), XML::xmlAttrs)[2]
  diseases <- data.frame(disease, disease.id)
  return(diseases)
}

#' Extract the Cellosaurus accession ID of one cell line.
#'
#' This function extract Cellosaurus accession ID of one cell line from an
#' XMLIntervalElementNode object extracted from Cellosaurus xml file. Combining
#' this function with \conde{\link[base]{apply}} or \code{\link[base]{sapply}}
#' can extract disease information from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A character. The Cellosaurus accession ID of cell line
#'
#' @examples
#' node <- GetAllCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
#' cell <- GetCell(node, "name", c("U251", "U87"))
#'
#' # get Cellosaurus Accession for first cell line
#' accession <- GetAccession(cell[[1]])
#'
#' # get Cellosaurus Accession for all cell lines
#' accessions <- sapply(cell, GetAccession)
GetAccession <- function(node){
  accession <- XML::xmlValue(XML::xmlChildren(node)$`accession-list`)
}

#' Check the Cellosaurus data version.
#'
#' Cellosaurus publishes its dataset as an XML file which could be downloaded
#' from ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml. This
#' could be used to check the publish version of this file (by default) or a
#' XML downloaded before (by passing the "file" argument).
#'
#' @param file A character. The xml file contains the Cellosaurus dataset.
#'
#' @examples
#' # Check the online Cellosaurus database version.
#' CellVersion()
#'
#' # Check the local Cellosaurus XML document version.
#' CellVersion(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
#'
#' @export
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
#' @param type A charatcer. It indicate the type of \code{id}. It could be
#' "names", "accession".
#'
#' @param id A vector of characters. It is the name or Cellosaurus accession of
#' cell lines that will be searched in Cellosaurus XML file.
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
GetCell <- function(node, type, id){
  if(type == "name"){
    xpath <- paste0("//name[text()='",
                    paste(id, collapse = "' or text() = '"),
                    "']/ancestor::cell-line")
  } else if (type == "accession") {
    xpath <- paste0("//accession[text()='",
                    paste(id, collapse = "' or text() = '"),
                    "']/ancestor::cell-line")
  }
    cells <- XML::getNodeSet(node, xpath)
}

#' Extract cell line information
#'
#' This function will extract the primary name, synonyms, Cellosausurs Accession
#' ID, disease, disease_id) of the cell lines and Wrap them into one data frame.
#' If you prefer some not all of these data, \code{\link{GetNames}},
#' \code{\link{GetDisease}}, or \code{\link{GetAccession}} are recommended.
#'
#' @param node An "XMLInternalElementNode" extracted from the Cellosaurus XML
#' file by either \code{\link{GetAllCell}} or \code{\link{GetCell}}
#'
#' @return A data frame contains 5 variables are:
#'   \item{name}{primary name of cell line.}
#'   \item{synonyms}{synonyms of cell lines separated by semicolons.}
#'   \item{accession}{Cellosaurus Accession ID for cell lines.}
#'   \item{disease}{diseases that are associated with the cell lines.}
#'   \item{disease_id}{NCI Thesaurus entry code of the disease.}
#'
#' @examples
#' all.cells <- GetAllCell(system.file("extdata",
#'                                     "cellosaurus.xml",
#'                                      package = "TidyComb"))
#' cell.lines <- GetCell(all.cells, c("U87", "A549"))
#' cell.info <- GetCellInfo(cell.lines)
#'
#' @export
GetCellInfo <- function(node) {
  stepi <- 1
  n <- length(node)

  names <- data.frame()
  diseases <- data.frame()
  accession <- character()
  temp <- NA

  for (i in 1:n) {
    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    flush.console()

    temp <- GetNames(node[[i]])
    names <- rbind.data.frame(names, temp)
    temp <- NA

    temp <- GetDisease(node[[i]])
    diseases <- rbind.data.frame(diseases, temp)
    temp <- NA

    temp <- GetAccession(node[[i]])
    accession <- c(accession, temp)
    temp <- NA

    stepi <- stepi + 1
  }
  info <- data.frame(names, accession, diseases)
}
