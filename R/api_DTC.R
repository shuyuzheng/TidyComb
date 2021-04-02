# TidyComb
# Functions for retrieving or updating drug information from  DTC database.

#' Match target according to CHEMBL IDs
#'
#' This function matches target of drugs according to user-provided identifiers.
#'
#' This function queries DTC database via
#' \href{chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/file:///Users/shuzhe/Desktop/Api_documentation.pdf}{API}
#' to query targets with CHEMBL ID
#'
#' @param chembl_ids A vector of characters contains the CHEMBL IDs of drugs
#'
#' @return A data frame contains two columns:
#' \itemize{
#'   \item \code{chembl_id} Input CHEMBL ID
#'   \item \code{target_gene_name} Gene name for the target
#'   \item \code{inhibitor_type} The type of the inhibitor
#'   \item \code{target_protein_class} The protein class for target
#'   \item \code{target_organism} The organism used to test the target activity
#'   \item \code{target_protein_name} The protein name for target
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' GetTarget("CHEMBL25")
GetTarget <- function(chembl_ids) {
  message("Getting targets from DTC...")
  url.base <- "https://drugtargetcommons.fimm.fi/api/data/bioactivity/?chembl_id="
  stepi <- 1
  n <- length(chembl_ids)
  target <- NULL
  for (id in chembl_ids) {
    temp <- NA
    tryCatch({

      url <- paste0(url.base, id)

      doc <- jsonlite::fromJSON(url)

      if (length(doc$bioactivities) != 0){
        temp <- doc$bioactivities
        temp$target_organism[is.na(temp$target_organism)] <- ""
        temp$target_pref_name[is.na(temp$target_pref_name)] <- ""
        temp$target_pref_name[temp$target_pref_name == "UNCHECKED"] <- ""
        temp <- temp[
          which(!tolower(temp$target_pref_name) == tolower(temp$target_organism)),
          c("chembl_id", "gene_name", "inhibitor_type", "protein_class",
            "target_organism", "target_pref_name")]
        colnames(temp) <- c(
          "chembl_id", "target_gene_name", "inhibitor_type",
          "target_protein_class", "target_organism", "target_protein_name")
        temp <- unique(temp)
        temp$target_organism[temp$target_organism == ""] <- NA
        temp$target_protein_name[temp$target_protein_name == ""] <- NA
      } else {
        temp <- data.frame(
          chembl_id = id,
          target_gene_name = NA,
          inhibitor_type = NA,
          target_protein_class = NA,
          target_organism = NA,
          target_protein_name = NA,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      temp <<- data.frame(
        chembl_id = id,
        target_gene_name = NA,
        inhibitor_type = NA,
        target_protein_class = NA,
        target_organism = NA,
        target_protein_name = NA,
        stringsAsFactors = FALSE
      )
    })
    target <- rbind.data.frame(target, temp)

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    stepi <- stepi + 1
  }
  return(target)
}
