# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

UnichemVer <- function() {
  url <- ""
}



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
#' @param ids A vector of characters. It contains identifiers of the drug.
#'
#' @param type A character or an integer. It specifies the type of identifiers
#' passed to \code{ids}. Acceptable values are:
#' \enumerate{
#'   \item \strong{An integer} indicates the "src_id" that UniChem assigned to
#'    various sources. A list of valid "src_id" can be found either on the
#'    \href{https://www.ebi.ac.uk/unichem/ucquery/listSources}
#'    {UniChem sources page} or by using \code{\link{GetAllSrcIds}} function.
#'   \item \strong{inchikey} Standard InChIKey of the drugs.
#' }
#'
#' @return A data frame with first column is the identifier passed to \code{ids}
#' argument and following columns are indentifiers of all resources integrated
#' by UniChem.
#'
#' @export
#'
#' @examples
#' aspirin <- GetIds(ids = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", type = "inchikey")
GetIds <- function(ids, type) {
  src.id <- GetAllSrcIds()
  # building base url
  if (type == "inchikey") {
    url.base <- "https://www.ebi.ac.uk/unichem/rest/inchikey/"
  } else if (type == "src.id")

  df <- data.frame() # result container
  i <- 1 # indicator
  n <- length(inchikey)

  # retrieving IDs

  for (id in inchikey) {
    # progress indicator
    message(round(i/n * 100), "%", "\r", appendLF = FALSE)
    flush.console()

    # requesting server via API
    url <- paste0(url.base, id) # establishing data
    res <- fromJSON(url) # parsing JSON file

    # extracting IDs
    drugbank <- ifelse(2 %in% res[,1], res[which(res[,1] == 2), 2], NA)
    kegg <- ifelse(6 %in% res[,1], res[which(res[,1] == 6), 2], NA)

    # assembling results
    temp <- data.frame(inchikey = id,
                       drugbank = drugbank,
                       kegg = kegg,
                       stringsAsFactors = FALSE)
    df <- rbind.data.frame(df, temp)

    # cleaning up temporary variable and updating indicator .
    temp <- data.frame()
    i <- i + 1
  }

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
