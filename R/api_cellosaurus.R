# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyright: Shuyu Zheng


#' Check Cellosaurus XML file version.
#'
#' \code{CellVersion} is used to check the version of a Cellosaurus XML file.
#'
#' Cellosaurus publishes its dataset as an XML file which could be downloaded
#' from \url{ftp://ftp.expasy.org/databases/cellosaurus/}. \code{CellVersion} is
#' used to check the dataset version of a "cellosaurus.xml" file. By default it
#' will check the online dataset version. If you'd like to check the version of
#' a local file, please pass the path of file to argument \code{file}.
#'
#' @param file File path to an XML file contains the Cellosaurus dataset.
#'
#' @param ... Other arguments required by \code{\link[XML]{xmlEventParse}}
#'
#' @return A named character vector. It contains: version, date of update,
#' number of archived cell lines and number of archived publications.
#'
#' @export
#'
#' @examples
#' # Check the online Cellosaurus dataset version.
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
  version <- NULL
  # Building handler
  startElement = function(ctxt, name, attrs, ...) {
    if (name == "release") {
      version <<- attrs
      XML::xmlStopParser(ctxt)
    }
  }
  class(startElement) = "XMLParserContextFunction"
  # Persing xml
  invisible(XML::xmlEventParse(file,
                               handlers = list(startElement = startElement)))
  return(version)
}

#' Update local Cellosaurus XML file
#'
#' \code{UpdateCell} checks and compares the versions of online and
#' local Cellosaurus XML file. If the local file is not up-to-date, it will
#' automatically update the local file.
#'
#' For parsing Cellosaurus XML file, we do \strong{not} recommend to directly
#' parse the online Cellosaurus dataset. It requires a huge amount of memory and
#' will lead to a crash of R. We recommend to download the dataset and then
#' parse the local XML file. In this case, it is necessary to update your local
#' file before processing it.
#'
#' \code{UpdateCell} using the \code{\link{CellVersion}} function to check the
#' current online version of Cellosaurus dataset and the version of a local
#' file. It will automatically update the local file if it is not up-to-date.
#'
#' Warning: The local file will be overwritted when updating data. If you'd like
#' to keep the old version of dataset, please backup it before updating.
#'
#' @param file The path to the local Cellosaurus file.
#'
#' @return A message about the version checking results. If the local file is
#' \strong{not} up-to-date, the local file will be updated with online
#' Cellosaurus data.
#'
#' @export
#'
#' @examples
#' UpdateCell(system.file("extdata", "cellosaurus.xml", package = "TidyComb"))
UpdateCell <- function(file) {
  version.local <- as.numeric(CellVersion(file)["version"])
  version.online <- as.numeric(CellVersion()["version"])
  if (version.local < version.online) {
    message("Online Cellosaurus dataset published new version ",
            version.online, "\n",
            " Local file is under version ", version.local, "\n",
            "Updating local files for you, please wait for a monment...")
    utils::download.file(url, file)
    message("Local Cellosaurus file hase been up to date now!")
  } else if (version.local == version.online) {
    message("Local file is already the latest version ",
            version.local, "No need for updating.")
  } else {
    message("There is something wrong with local Cellosaurus file. ",
            "Please check it.")
  }
}

#' Loading a Cellosaurus XML dataset
#'
#' \code{GetAllCell} will parse the Cellosaurus XML file and extract all content
#' in "cell-line-list" node as a \code{XML document} object.
#'
#' Cellosaurus.xml file contains 5 child nodes in its root node:
#' "header", "cell-line-list", "publication-list", "copyright". (more
#' information in "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xsd")
#' All the cell line informations we need for preparing data are in
#' "cell-line-list" so this function will parse the dataset file and remove all
#' rudundant informations.
#'
#' Warning: Although it is possible to parsing the online database directly by
#' passing url \url{ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml}
#' to function \code{GetAllCell}, it is easy to cause crach of R as the huge
#' reuirement of memory. We recommend to download the dataset to a local file
#' and then parse this local file by using this function.
#'
#' @param file File path to a Cellosaurus xml file.
#'
#' @return An XMLNode containing all cell lines' information archieved in
#' Cellosaurus dataset.
#'
#' @export
#'
#' @examples
#' all.cell <- GetAllCell(system.file("extdata", "cellosaurus.xml",
#'                                    package = "TidyComb"))
GetAllCell <- function(file) {
  doc <- XML::xmlInternalTreeParse(file)
  all.cell <- XML::xmlRoot(doc)[[2]]
  return(all.cell)
}

#' Find matching cell-lines
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
#' @export
#'
#' @examples
#' node <- GetAllCell(file = system.file("extdata", "cellosaurus.xml",
#'                                       package = "TidyComb"))
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
#' \code{GetCellName} extract primary name and synonyms from only one
#' \emph{cell-line-list} node in Cellosaurus xml file.
#'
#' This function extracts the primary name (as "name") and synonyms (as
#' "synonyms") of an \code{XMLIntervalElementNode} object containing
#' one \emph{cell-line-list} node extracted from Cellosaurus xml file.
#'
#' If you'd like to extract information from multiple \emph{cell-line-list}
#' nodes, combining this function with \code{xpathSapply} or
#' \code{sapply} is recommanded.
#'
#' @param node An \code{XMLInternalElementNode} with only one
#' \emph{cell-line-list} node extracted from \emph{Cellosaurus xml file}.
#'
#' @return A data frame contains two variables:
#' \itemize{
#'   \item \code{name} the primary name of cell line.
#'   \item \code{synonyms} synonyms of the cell line separated by semicolons.
#' }

GetCellName <- function(node) {
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

#' Extract primary name and synonyms of one cell line into one character.
#'
#' \code{GetCellNameInOne} extract primary name and synonyms from only one
#' \emph{cell-line-list} node in Cellosaurus xml file and paste them into one
#' character.
#'
#' This function extracts all names(primary name and synonyms) of an
#' \code{XMLIntervalElementNode} object containing one \emph{cell-line-list}
#' node extracted from Cellosaurus xml file. Then paste them into one character.
#'
#' If you'd like to extract information from multiple \emph{cell-line-list}
#' nodes, combining this function with \code{xpathSapply} or \code{sapply} is
#' recommanded.
#'
#' @param node An \code{XMLInternalElementNode} with only one
#' \emph{cell-line-list} node extracted from \emph{Cellosaurus xml file}.
#'
#' @return A character contains all names (Primary and synonyms) of cell line.

GetCellNameInOne <- function(node) {
  name.list <- XML::xmlToDataFrame(XML::xmlChildren(node)$'name-list',
                                   stringsAsFactors = FALSE)
  names <- sapply(name.list, paste, collapse = "; ")
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
#' this function with \code{apply} or \code{sapply} can extract disease
#' information from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A data frame contains two variables:
#' \item{disease}{name of the cell line associated disease.}
#' \item{disease_id}{NCI Thesaurus entry code of the associated disease.}
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
#' Combining it with \code{apply} or \code{sapply} can extract tissue from an
#' XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A character. The tissue name of cell line according to CClE category.
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
#' this function with \code{apply} or \code{sapply} can extract disease
#' information from an XMLNodeSet object.
#'
#' @param node An XMLInternalElementNode with only one cell line's information
#' which was extracted from Cellosaurus xml file.
#'
#' @return A character. The Cellosaurus accession ID of cell line
GetAccession <- function(node){
  accession <- XML::xpathSApply(node,
                                './accession-list/accession[@type="primary"]',
                                XML::xmlValue)
}

#' Extract cell line information
#'
#' This function will extract the primary name, synonyms, Cellosausurs Accession
#' ID, disease, disease_id) of the cell lines and Wrap them into one data frame.
#' If you prefer some not all of these data, \code{\link{GetCellName}},
#' \code{\link{GetCellNameInOne}}, \code{\link{GetDisease}}, or
#' \code{\link{GetAccession}} are recommended.
#'
#' @param node An "XMLInternalElementNode" extracted from the Cellosaurus XML
#' file by either \code{\link{GetAllCell}} or \code{\link{GetCell}}
#'
#' @param info A Character indicate about what kind of infomation you want to
#' extract. Available values are:
#' \itemize{
#'   \item \strong{name} primary name and synonyms of cell lines.
#'   \item \strong{name_in_one} Output all names (Primary and synonyms) of cell
#'   line into one character.
#'   \item \strong{accession} Cellosaurus Accession ID for cell lines.
#'   \item \strong{disease} diseases that are associated with the cell lines and
#'    corresponding NCI Thesaurus entry code.
#'   \item \strong{tissue} Tissues that the cell line generated from. Following
#'   the CCLE category.
#'}
#'
#' @return A data frame contains cell line information selected by \code{info}
#'
#' @export
#'
#' @examples
#' node <- GetAllCell(system.file("extdata", "cellosaurus.xml",
#'                                package = "TidyComb"))
#' cell.lines <- GetCell(node, c("U87", "A549"), "name")
#' cell.info <- GetCellInfo(cell.lines)
GetCellInfo <- function(node, info = "accession") {

  if (info == "name") {
    fun <- function(x) {GetCellName(x)}
    colname <- c("name", "synonyms")
  } else if (info == "name_in_one") {
    fun <- function(x) {GetCellNameInOne(x)}
    colname <- c("all_names")
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
    stop("Info ", info, ' is not allowed. Available values are: "name",',
         '"accession", "disease", "tissue", and "name_in_one".')
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

