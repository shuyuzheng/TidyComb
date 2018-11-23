library(XML)
library(RCurl)

GetChembl <- function(ids, quiet = TRUE){
  # Input: character vector of compound names
  # Output: Data frame with matched Chembl ID, max clinical phase.
  #
  # API Documentation:https://www.ebi.ac.uk/chembl/ws

  curlHandle <- getCurlHandle()
  out <- data.frame(stringsAsFactors = FALSE)
  new_item <- NA
  
  stepi <- 0
  n <- length(ids)
  for (compound in ids) {
    
    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    flush.console()
    
    tryCatch({
      chembl_phase <- integer()
      chembl_id <- character()
      res <- dynCurlReader()
      url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", compound)
      curlPerform(url = url, curl = curlHandle, writefunction = res$update)
      doc <- xmlInternalTreeParse(res$value())
      new_item <- xpathApply(doc, "//max_phase", xmlValue)
      new_item <- as.integer(unlist(new_item))
      if (is.null(new_item)) {
        new_item <- NA
      }
      chembl_phase <- c(chembl_phase, new_item)
     
      new_item <- xpathApply(doc, "//molecule_hierarchy/molecule_chembl_id", xmlValue)
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
      chembl_id <- c(chembl_id, new_item)

      },
    error = function(e) {
      chembl_id <<- NA
      chembl_phase <<- 0
      print(e)
    }
    )
    df <- data.frame(id = compound, chembl_id = chembl_id,
                     chembl_phase = chembl_phase,
                     stringsAsFactors = F)
    out <- rbind.data.frame(out, df)
    stepi <- stepi + 1
  }

    # Cleanup
  #rm(curlHandle)
  #gc()
  return(out)
}
