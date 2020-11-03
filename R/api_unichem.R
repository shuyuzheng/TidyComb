################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng<shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for retrieving or updating drug information from UniChem.

#' Get other identifiers from InChIKey
#'
#' \code{GetIds} matches the identifiers of
#' \href{https://www.ebi.ac.uk/unichem/ucquery/listSources}{UniChem sources},
#' including DrugBank, KEGG, Zinc...
#'
#' This function retrieves the of drug
#' via \href{https://www.ebi.ac.uk/unichem/info/webservices}{UniChem Web
#' Services}.
#'
#' @param inchikey A vector of characters. It contains identifiers of the drug.
#'
#' @return A data frame with first column is the identifier passed to \code{ids}
#' argument and following columns are indentifiers of all resources integrated
#' by UniChem.
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' aspirin <- GetIds(inchikey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
GetIds <- function(inchikey) {
  message("Getting identifiers from UniChem...")
  # building base url
    url.base <- "https://www.ebi.ac.uk/unichem/rest/inchikey/"

  df <- NULL
  i <- 1 # indicator
  n <- length(inchikey)

  # retrieving IDs

  for (id in inchikey) {
    # progress indicator
    message(round(i/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    drugbank <- NA
    kegg <- NA
    chembl <- NA
    tryCatch({
      # requesting server via API
      url <- paste0(url.base, id) # establishing data
      res <- jsonlite::fromJSON(url) # parsing JSON file

      # extracting IDs
      if (2 %in% res[, 1]) {
        drugbank <- res[which(res[, 1] == 2), 2]
      }
      if (6 %in% res[, 1]) {
        kegg <- res[which(res[, 1] == 6), 2]
      }
      if (1 %in% res[, 1]) {
        chembl <- res[which(res[, 1] == 1), 2]
      }

    }, error = function(e) {
      print(e)
    }
    )
    temp <- cbind(inchikey = id,
                  uni_drugbank = drugbank,
                  uni_kegg_c = kegg,
                  chembl_id = chembl)
    df <- rbind(df, temp)
    i <- i + 1
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  return(df)
}


#' Get all src_ids
#'
#' Obtain all resources and corresponding \code{src_ids} currently in UniChem
#'
#' @return A data frame contains 3 columns:
#' \item{name}{Name of the resources.}
#' \item{src_id}{"src_id" of the resources.}
#' \item{description}{A brief description of the content of the source.}
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' resource <- GetAllSrcIds()
GetAllSrcIds <- function(){
  url <- "https://www.ebi.ac.uk/unichem/rest/src_ids/"
  src.ids <- jsonlite::fromJSON(url)[,1]
  df <- data.frame()

  url.base <- "https://www.ebi.ac.uk/unichem/rest/sources/"
  for (id in src.ids) {
    temp <- data.frame()

    url <- paste0(url.base, id)
    temp <- jsonlite::fromJSON(url)[c("name", "src_id", "description")]
    df <- rbind.data.frame(df, temp)
  }
  return(df)
}
