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
#'
#' @param ... Other arguments.
#'
#' @export
#'
#' @examples
#' # Check the online Cellosaurus database version.
#' CellVersion()
#'
#' # Check the local Cellosaurus XML document version.
#' CellVersion(system.file("extdata",
#'                         "cellosaurus.xml",
#'                         package = "TidyComb"))
#'

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
#' @return An XMLNode containing all cell lines' information archieved in
#' Cellosaurus dataset.
#'
#' @export
#'
#' @examples
#' all.cell <- GetAllCell()
GetAllCell <- function() {
  doc <- XML::xmlInternalTreeParse(system.file("extdata",
                                               "cellosaurus.xml",
                                               package = "TidyComb"))
  all.cell <- XML::xmlRoot(doc)[[2]]
  return(all.cell)
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
#' @param ids A vector of characters. It is the name or Cellosaurus accession of
#' cell lines that will be searched in Cellosaurus XML file.
#'
#' @param type A charatcer. It indicate the type of \code{id}. It could be
#' "name", "accession".
#'
#' @return An XMLNodeSet containing matched cell lines in the dataset. If no
#' cell line is matched, a NULL list will be return.
#'
#' @examples
#' node <- GetAllCell()
#' cell.lines <- GetCell(node, c("U87", "A549"), "name")
#'
GetCell <- function(node, ids, type = "name"){
  if(type == "name"){

    xpath <- paste0("//name[text()='",
                    paste(ids, collapse = "' or text() = '"),
                    "']/ancestor::cell-line")
  } else if (type == "accession") {
    xpath <- paste0("//accession[text()='",
                    paste(ids, collapse = "' or text() = '"),
                    "']/ancestor::cell-line")
  } else {
    stop("Type ", type,
         'is not allowed. Available types are: "name", "accession"')
  }
    cells <- XML::getNodeSet(node, xpath)
}

#' Extract primary name and synonyms of one cell line.
#'
#' \code{GetNames} extract primary name and synonyms from only one
#' \emph{cell-line-list} node in Cellosaurus xml file.
#'
#' This function extracts the primary name (as "name") and synonyms (as
#' "synonyms") of an \code{XMLIntervalElementNode} object containing
#' one \emph{cell-line-list} node extracted from Cellosaurus xml file.
#'
#' If you'd like to extract information from multiple \emph{cell-line-list}
#' nodes, combining this function with \code{xpathSapply} or
#' \code{\link[base]{sapply}} is recommanded.
#'
#' @param node An \code{XMLInternalElementNode} with only one
#' \emph{cell-line-list} node extracted from \emph{Cellosaurus xml file}.
#'
#' @return A data frame contains two variables:
#' \itemize{
#'   \item \code{name} the primary name of cell line.
#'   \item \code{synonyms} synonyms of the cell line separated by semicolons.
#' }
#'
#' @examples
#' # parse the Cellosaurus xml file.
#' node <- GetAllCell()
#'
#' # extract "cell-line-list" nodes of "U251" and "U87" cell lines.
#' cell <- GetCell(node, c("U251", "U87"), "name")
#'
#' # get first cell line (U251)'s primary name and synonyms
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
  names <- cbind(name = name, synonyms = synonyms)
  return(names)
}

#' Extract the disease information of one cell line.
#'
#' \code{GetDisease} extracts the cell line associated diseases and NCI
#' Thesaurus entry code from one \code{cell-line-list} node in Cellosaurus xml
#' file.
#'
#' This function extract the cell line associated disease and disease ID (NCI
#' Thesaurus entry code) which the cell line associated with from an
#' XMLIntervalElementNode object extracted from Cellosaurus xml file. Combining
#' this function with \code{\link[base]{apply}} or \code{\link[base]{sapply}}
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
#' node <- GetAllCell()
#' cell <- GetCell(node, c("U251", "U87"), "name")
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
  diseases <- cbind(disease, disease.id)
  return(diseases)
}

#' Extract the source tissue of one cell line.
#'
#' This function extract source tissue according to cross-reference "CCLE Name"
#' from an XMLIntervalElementNode object extracted from Cellosaurus xml file.
#' Combining it with \code{\link[base]{apply}} or \code{\link[base]{sapply}}
#' can extract tissue from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A character. The tissue name of cell line according to CClE category.
#'
#' @examples
#' node <- GetAllCell()
#' cell <- GetCell(node, c("U251", "U87"), "name")
#'
#' # get Cellosaurus Accession for first cell line
#' accession <- GetAccession(cell[[1]])
#'
#' # get Cellosaurus Accession for all cell lines
#' accessions <- sapply(cell, GetAccession)
GetTissue <- function(node){
  ref.list <- XML::xmlChildren(node)$`xref-list`
  ref <- sapply(XML::xmlChildren(ref.list), XML::xmlAttrs)
  ccle <- ref[3, which(ref[1,] == "CCLE")]
  tissue <- tolower(gsub("^[^_]+(?=_)_", "",ccle, perl = TRUE))
}

#' Extract the Cellosaurus accession ID of one cell line.
#'
#' This function extract Cellosaurus accession ID of one cell line from an
#' XMLIntervalElementNode object extracted from Cellosaurus xml file. Combining
#' this function with \code{\link[base]{apply}} or \code{\link[base]{sapply}}
#' can extract disease information from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A character. The Cellosaurus accession ID of cell line
#'

GetAccession <- function(node){
  accession <- XML::xpathSApply(node,
                                './accession-list/accession[@type="primary"]',
                                XML::xmlValue)
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
#' @param info A Character indicate about what kind of infomation you want to
#' extract. Available values are:
#' \itemize{
#'   \item \strong{name} primary name and synonyms of cell lines.
#'   \item \strong{accession} Cellosaurus Accession ID for cell lines.
#'   \item \strong{disease} diseases that are associated with the cell lines and
#'    corresponding NCI Thesaurus entry code.
#'   \item \strong{tissue} Tissues that the cell line generated from. Following the CCLE
#'    category.
#'}
#'
#' @return A data frame contains cell line information selected by \code{info}
#'
#' @examples
#' node <- GetAllCell()
#' cell.lines <- GetCell(node, c("U87", "A549"), "name")
#' cell.info <- GetCellInfo(cell.lines)
#'
#' @export
GetCellInfo <- function(node, info = "accession") {

  if (info == "name") {
    fun <- function(x) {GetNames(x)}
    colname <- c("name", "synonyms")
  } else if (info == "accession") {
    fun <- function(x) {GetAccession(x)}
    colname <- c("cellosaurus_accession")
  } else if (info == "disease") {
    fun <- function(x) {GetDisease(x)}
    colname <- c("disease_name", "disease_id")
  } else if (info == "tissue") {
    fun <- function(x) {GetTissue(x)}
    colname <- c("tissue_name")
  } else {
    stop("Info ", info, 'is not allowed. Available values are: "name",',
         '"accession", "disease", "tissue".')
  }

  temp <- NULL
  mat <- NULL

  stepi <- 1
  n <- length(node)
  for (i in 1:n) {
    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    temp <- fun(node[[i]])
    if (length(temp) == 0) {
      temp <- rep(NA, length(colname))
    }
    mat <- rbind(mat, temp)
    temp <- NULL

    stepi <- stepi + 1
  }

  colnames(mat) <- colname
  return(mat)
}

