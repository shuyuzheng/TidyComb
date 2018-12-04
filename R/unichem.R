# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

UnichemVer <- function() {
  url <- ""
}



#' Extract information from Chembl database
#'
#' This function retrieves the Chembl ID and max clinical trial phase of drug
#' via \href{https://www.ebi.ac.uk/unichem/info/webservices}{UniChem Web
#' Services}.
#'
#' @param ids A vector of characters. It contains identifiers of the drug.
#'
#' @param type A character. It specifies the type of identifiers passed to
#' \code{ids}. Acceptable inputs are:
#' \enumerate{
#'   \item src_compound_id: It is the individual compound identifier provided by
#' each of the sources.
#'   \item src_id: It is integers that represent the various sources in UniChem.
#' A list of valid src_id's can be found either on the
#' \href{https://www.ebi.ac.uk/unichem/ucquery/listSources}
#' {UniChem sources page} or by using \code{\link{GetAllSrcIds}} function.
#'   \item InChIKey:
#' }
#'
#' @return
#'
GetUnichem <- function(ids, type) {

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
    temp <- fromJSON(url)[c("name", "src_id", "description")]
    df <- rbind.data.frame(df, temp)
  }
  return(df)
}
