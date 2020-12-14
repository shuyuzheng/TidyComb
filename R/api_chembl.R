################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for retrieving or updating drug information from ChEMBL database.
#
# Functions on this page
# ChembelVersion: Check the version of ChEMBL databse
# GetChembl: Get informateion of drugs from ChEMBL
# GetChemblPhase: Get clinical trial phase of drugs from ChEMBL

#' Check version of ChEMBL database
#'
#' By calling function \code{ChemblVersion} without any argument, you can check
#' the current version ChEMBL database.
#'
#' \code{ChemblVersion} is a fuction wrapping the "Status" function provided by
#' \href{https://www.ebi.ac.uk/chembl/ws}{ChEMBL Web Services}.
#'
#' @return A vector of named characters contains the update information of
#' ChEMBL.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' # print(ChemblVersion())
ChemblVersion <- function(){
  url <- "https://www.ebi.ac.uk/chembl/api/data/status?format=json"
  res <- jsonlite::fromJSON(url, flatten = TRUE)
  res <- unlist(res)
  return(res)
}

#' Get information from ChEMBL
#'
#' \code{GetChembl} retrieves ChEMBL ID and max clinical trial phase of drugs
#' accordint to standard InChIKey.
#'
#' \code{GetChembl} queries ChEMBL database via
#' \href{https://www.ebi.ac.uk/chembl/ws}{ChEMBL Web Services}. The input drug
#' identifier is standard InChIKey.
#'
#' @param ids A vector of characters contains the InChIKey of drugs.
#' \emph{Note:} the leading and tailing whitespaces are not allowed.
#'
#' @param info A vector of characters contains the information user want to get
#' from ChEMBL database. The available values are:
#'  \item{"phase"} {Maximun clinical phase for input compounds.}
#'  \item{"molecular_formula"} {Molecula formula for input compounds.}
#'  \item{"inchikey"} {Standard InChIKey for input compounds.}
#'  \item{"smiles"} {Canonical SMILES for input compounds.}
#'  \item{"chembl_id"} {ChEMBL ID for input compounds.}
#'
#' @param quiet A logical value. If it is \code{TRUE}, the error messages will
#' not show during runing the function.
#'
#' @return A data frame contains reqired information for drugs and one column
#' "input_id" for input identifier.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' drug.info <- GetChembl("PMATZTZNYRCHOR-CGLBZJNRSA-N",
#'                        info = c("phase", "chembl_id", "inchikey"))
GetChembl <- function(ids, info = c("chembl_id", "phase"), quiet = TRUE) {
  message("Getting information from ChEMBL...")
  vail_info <- c("phase", "molecular_formula", "chembl_id", "inchikey", "smiles")
  invail_info <- setdiff(info, vail_info)
  if (length(invail_info) > 0){
    stop("The required info '", paste0(invail_info, collapse = ", "), "' is/are",
    "not vailable for this function. A vailable values are")
  }
  curlHandle <- RCurl::getCurlHandle()
  out <- NULL

  stepi <- 1
  n <- length(ids)
  for (id in ids) {

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    tryCatch({
      df <- NULL
      res <- RCurl::dynCurlReader()

      url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", id)
      RCurl::curlPerform( url = url,
                          curl = curlHandle, writefunction = res$update)
      doc <- XML::xmlInternalTreeParse(res$value())

      item <- sapply(info, function(x){
        field_path <- switch(x,
          phase = "//molecule/max_phase",
          molecular_formula = "//molecule/molecule_properties/full_molformula",
          chembl_id = "//molecule/molecule_chembl_id",
          inchikey = "//molecule/molecule_structures/standard_inchi_key",
          smiles = "//molecule/molecule_structures/canonical_smiles",
        )
        return(.extinfoXML(doc, field_path))
      })
      # clinical phase
      # chembl_phase <- XML::xpathApply(doc, "//max_phase", XML::xmlValue)
      # chembl_phase <- as.integer(unlist(chembl_phase))
      # if (is.na(chembl_phase)) {
      #   chembl_phase <- 0
      # }
      #
      # # Chembl ID
      # chembl_id <- XML::xpathApply(doc, "//molecule/molecule_chembl_id",
      #                             XML::xmlValue)
      # chembl_id <- unlist(chembl_id)
      # if (is.null(chembl_id)) {
      #   chembl_id <- NA
      # }
      #
      # # InChIKey
      # inchikey <- XML::xpathApply(doc, "//molecule/molecule_structures/standard_inchi_key",
      #                              XML::xmlValue)
      # inchikey <- unlist(inchikey)
      # if (is.null(inchikey)) {
      #   inchikey <- NA
      # }
      #
      # # smiles
      # inchikey <- XML::xpathApply(doc, "//molecule/molecule_structures/standard_inchi_key",
      #                             XML::xmlValue)
      # inchikey <- unlist(inchikey)
      # if (is.null(inchikey)) {
      #   inchikey <- NA
      # }
      # # InChIKey
      # inchikey <- XML::xpathApply(doc, "//molecule/molecule_structures/standard_inchi_key",
      #                             XML::xmlValue)
      # inchikey <- unlist(inchikey)
      # if (is.null(inchikey)) {
      #   inchikey <- NA
      # }
    },
    error = function(e) {
      item <<- sapply(info, function(x){
        return(assign(x, NA))
      })
      if (!quiet) {
        print(e)
      }
    })
    df <- as.data.frame(as.list(item))
    df$input_id <- id
    # df <- data.frame(inchikey = id, chembl_id = chembl_id,
    #                  chembl_phase = chembl_phase,
    #                  stringsAsFactors = FALSE)
    out <- rbind.data.frame(out, df)
    stepi <- stepi + 1
  }

  # Cleanup
  rm(curlHandle)
  gc()

  return(out)
}

#' Get max clinical trial phase from ChEMBL ID
#'
#' \code{GetChemblPhase} retrieves max clinical trial phase of drugs accordint
#' to ChEMBL ID.
#'
#' \code{GetChemblPhase} queries ChEMBL database via
#' \href{https://www.ebi.ac.uk/chembl/ws}{ChEMBL Web Services}. The input drug
#' identifier is ChEMBL ID and the out put information is "max clinical trial
#' phase".
#'
#' @param ids A vector of characters contains the ChEMBL IDs of drugs.
#' \emph{Note:} the leading and tailing whitespaces are not allowed.
#'
#' @param quiet A logical value. If it is \code{TRUE}, the error messages will
#' not show during runing the function.
#'
#' @return A data frame contains 2 columns:
#'  \item{chembl_id}{Inputted ChEMBL ID of drugs.}
#'  \item{chembl_phase}{Max clinical trial phase of the matched drugs.}
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' drug.info <- GetChemblPhase("CHEMBL25")

GetChemblPhase <- function(ids, quiet = TRUE) {
  message("Getting clinical phase from ChEMBL...")
  chembl_id <- NULL
  chembl_phase <- NULL
  curlHandle <- RCurl::getCurlHandle()

  stepi <- 1
  n <- length(ids)
  for (id in ids) {

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    tryCatch({
      new_item <- 0
      res <- RCurl::dynCurlReader()

      url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", id)
      RCurl::curlPerform( url = url,
                          curl = curlHandle, writefunction = res$update)
      doc <- XML::xmlInternalTreeParse(res$value())

      # clinical phase
      new_item <- XML::xpathApply(doc, "//max_phase", XML::xmlValue)
      new_item <- as.integer(unlist(new_item))
      if (is.na(new_item)) {
        new_item <- 0
      }

    }, error = function(e) {
      if (!quiet) {
        new_item <<- 0
        print(e)
      }
    }
    )
    chembl_phase <- c(chembl_phase, new_item)

    chembl_id <- c(chembl_id, id)
    stepi <- stepi + 1
  }

  df <- data.frame(chembl_id, chembl_phase, stringsAsFactors = F)

  # Cleanup
  rm(curlHandle)
  gc()

  return(df)
}

.extinfoXML <- function(doc, field_path){
  new_item <- XML::xpathApply(doc, field_path,
                              XML::xmlValue)
  new_item <- unlist(new_item)
  if (is.null(new_item)) {
    new_item <- NA
  } else if (length(new_item) > 1){
    new_item <- paste0(new_item, collapse = "; ")
  }

  return(new_item)
}