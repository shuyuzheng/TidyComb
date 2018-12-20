# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

#' Check version of Chembl database
#'
#' @return A vector of named characters contains the update information of
#' Chembl.
#' @export
#'
#' @examples
#' print(ChemblVersion())
ChemblVersion <- function(){
  url <- "https://www.ebi.ac.uk/chembl/api/data/status?format=json"
  res <- jsonlite::fromJSON(url, flatten = TRUE)
  res <- unlist(res)
}

#' Get Chembl ID and max clinical phase
#'
#' This function retrieves the Chembl ID and max clinical trial phase of drugs
#' by provided InChiKeys.
#' API: \href{https://www.ebi.ac.uk/chembl/ws}{ChEMBL Web Services}.
#'
#' @param ids A vector of characters contains the InChIKey of drugs. Note: the
#' leading and tailing whitespaces are not allowed.
#'
#' @param quiet A logical value. If it is \code{TRUE}, the error messages will
#' not show during runing the function.
#'
#' @return A data frame contains 3 columns:
#'  \item{inchikey}{Identifiers inputted in the function.}
#'  \item{chembl_id}{Chembl ID of the matched drugs.}
#'  \item{chembl_phase}{Max clinical trial phase of the matched drugs.}
#'
#' @export
#'
#' @examples
#' drug.info <- GetChembl("PMATZTZNYRCHOR-CGLBZJNRSA-N")
GetChembl <- function(ids, quiet = TRUE){
  message("Getting information from ChEMBL...")
  curlHandle <- RCurl::getCurlHandle()
  out <- data.frame(stringsAsFactors = FALSE)
  new_item <- NA

  stepi <- 1
  n <- length(ids)
  for (id in ids) {

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    tryCatch({

      chembl_phase <- integer()
      chembl_id <- character()
      res <- RCurl::dynCurlReader()

      url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", id)
      RCurl::curlPerform(
        url = url,
        curl = curlHandle, writefunction = res$update)
      doc <- XML::xmlInternalTreeParse(res$value())

      # clinical phase
      new_item <- XML::xpathApply(doc, "//max_phase", XML::xmlValue)
      new_item <- as.integer(unlist(new_item))
      if (is.null(new_item)) {
        new_item <- 0
      }
      chembl_phase <- c(chembl_phase, new_item)

      # Chembl ID
      new_item <- XML::xpathApply(doc, "//molecule/molecule_chembl_id",
                                  XML::xmlValue)
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
      chembl_id <- c(chembl_id, new_item)
      },
    error = function(e) {
      chembl_id <<- NA
      chembl_phase <<- 0
      if (!quiet) {
        print(e)
      }
    }
    )
    df <- data.frame(inchikey = id, chembl_id = chembl_id,
                     chembl_phase = chembl_phase,
                     stringsAsFactors = F)
    out <- rbind.data.frame(out, df)
    stepi <- stepi + 1
  }

  # Cleanup
  rm(curlHandle)
  gc()

  return(out)
}