################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng<shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for retrieving or updating publication information from PubMed.

#'Get publication information from Pubmed
#'
#' \code{GetPubmed} extract information (pulication date, first name of first
#' author, title of the publication, full name of journal, doi) of publications
#' from Pubmed database.
#'
#' \code{GetPubmed} extract publications information via API provided by PMC
#' database.
#'
#' The argument \code{tool} and \code{email} are required for access
#' \href{https://www.ncbi.nlm.nih.gov/pmc/tools/developers/}{API}.
#' Following is the "Note about API usage":
#'
#' The use of our APIs is entirely free, and doesn't require an API key, but
#' we ask that you please:
#' \itemize{
#'   \item{Do not make concurrent requests, even at off-peak times;}
#'   \item{Include two parameters that help to identify your service or
#' application to our servers:
#'     \enumerate{
#'       \item tool should be the name of the application, as a string value
#'       with no internal spaces.
#'       \item email should be the e-mail address of the maintainer of the tool,
#'       and should be a valid e-mail address.
#'     }
#'   }
#' }
#'
#' @param pmids A vector of character/numeric containing PubMed ID of
#' publications
#'
#' @param tool The name of application, default is "TidyComb"
#'
#' @param email User's e-mail address
#'
#' @return A data frame containing 6 variables:
#' \itemize{
#'   \item \strong{pubmed_id} The PubMed ID of the publication passed to
#'   \code{pmid} argument.
#'   \item \strong{year} The year of the publication date.
#'   \item \strong{name} The first name of the first author.
#'   \item \strong{title} The title of the publication.
#'   \item \strong{journal} The full name of journal where the publication was
#'   reposted.
#'   \item \strong{doi} The Digital Object Identifier (DOI) of the publication.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' study <- GetPubmed(c(26551875, 15627163), email = "shuyu.zheng1992@gmail.com")
GetPubmed <- function(pmids, tool = "TidyComb", email = NULL ) {
  if (is.null(email)) {
    stop("As the PMC API ask for user's contact information, please provide
         your email address as a parameter when calling this function.")
  }

  df <- NULL
  for (id in pmids){
    url  <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                   "esummary.fcgi?", "db=pubmed", "&id=", id, "&tool=R",
                   "&email=", email, "&retmode=json")
    res <- jsonlite::fromJSON(url)
    result <- unlist(res[[2]], recursive = TRUE, use.names = TRUE)
    info <- as.data.frame(t(result[paste0(id, ".", c("pubdate", "authors.name1",
                                                     "title", "fulljournalname",
                                                     "elocationid"))]))
    colnames(info) <- c("year", "name", "title", "journal", "doi")
    info$year <- as.integer(substr(info$year, 1, 4))
    info$name <- sub("([A-Za-z]+).*", "\\1", info$name)
    info$doi <- sub("^.*doi: ", "", info$doi)
    info$pubmed_id <- as.integer(id)
    df <- rbind.data.frame(df, info)
  }

  return(df)
}

