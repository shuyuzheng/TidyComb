# TidyComb
# Functions for retrieving or updating celline information from Cellosaurus.
# Copyrighte Shuyu Zheng

#' Match CID according to other identifiers
#'
#' \code{GetCid} matches CIDs of drugs according to user-provided identifiers.
#'
#' This function using the PubChem API
#' (\href{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}{PUG REST}
#' to seach mached drugs. Available
#' identifiers are: name, SID, SMILES, InChI, SDF, InChIKey, molecula formula.
#'
#' \strong{API USAGE POLICY:}
#' Please note that PUG REST is not designed for very large volumes (millions)
#' of requests. We ask that any script or application not make more than 5
#' requests per second, in order to avoid overloading the PubChem servers. If
#' you have a large data set that you need to compute with, please contact us
#' for help on optimizing your task, as there are likely more efficient ways to
#' approach such bulk queries.
#'
#' @param ids A vector of characters contains the identifiers of drugs you are
#' interested in.
#'
#' @param type A character indicates the type of identifiers passed by \code{id}
#'   argument. Available types are:
#'   \itemize{
#'   \item \strong{smiles} Identifiers of drugs in the simplified molecular-input
#'   line-entry system (SMILES).
#'   \item \strong{inchi} The IUPAC International Chemical Identifier (InChI).
#'   \item \strong{inchikey} Standard InChIKey of the drugs.
#'   \item \strong{formula} Molecule formula of the drugs.
#'   \item \strong{name} Name for drugs. This type could be used for searching
#'   synonyms, NCGC IDs, Chembl IDs, CAS or any other kind of identifiers.
#'   \item \strong{sdf} The SDF
#'   }
#' @param quiet
#'
#' @return A data frame contains two columns:
#' \itemize{
#'   \item \code{id} The identifiers inputted by user.
#'   \item \code{cid} The matched CID.
#' }
#'
#' @export
#'
#' @examples
#' GetCid("aspirin", "name", quiet = TRUE)
GetCid <- function(ids, type = NULL, quiet = TRUE){

  types <- c("name", "smiles", "inchi", "sdf", "inchikey", "formula")

  if (!(type %in% types)) {
    stop("Invalid idtype specified, valiable idtypes are: ", types)
  }

  curlHandle <- getCurlHandle()
  out <- data.frame(stringsAsFactors = FALSE)

  stepi <- 0
  n <- length(ids)
  for (id in ids) {

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    flush.console()

    tryCatch({
      res <- dynCurlReader()
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound",
                    "/%s/%s", "/synonyms/XML")
      curlPerform(
        url = sprintf(url, type, URLencode(id)),
        curl = curlHandle, writefunction = res$update)
      doc <- xmlInternalTreeParse(res$value())
      rootNode <- xmlName(xmlRoot(doc))
      if (rootNode == "InformationList") {
        xpathApply(doc, "//x:Information", namespaces = "x", function(x) {
          cid <- xpathSApply(x, "./x:CID", namespaces = "x", xmlValue)
          df <- data.frame(id = id, cid = cid,
                           stringsAsFactors = FALSE)
          out <<- rbind.data.frame(out, df)
        })
      } else if (rootNode == "Fault") {
        xpathApply(doc, "//x:Details", namespaces = "x", function(x){
        fault <- xpathApply(doc, "//x:Details", namespaces = "x", xmlValue)
        if (!quiet) {
          print( paste(id, fault[[1]], sep = ": ") )
        }
        df <- data.frame(id = id, cid = NA,
                         stringsAsFactors = FALSE)
        out <<- rbind.data.frame(out, df)
       })
      }
    },
    error = function(e) {
      print(e)
    },
    finally = Sys.sleep(0.2) # See usage policy.
    )
    stepi <- stepi + 1
  }

    # Cleanup
  rm(curlHandle)
  gc()
  return(out)
}

#' Title
#'
#' @param cids
#'
#' @return
#' @export
#'
#' @examples
GetPubNames <- function(cids){

  url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     "%s", "/synonyms/JSON")

  # Build containers
  cid <- NULL
  name <- NULL
  synonyms <- NULL
  df <- data.frame(cid = integer(), name = character(), synonyms = character())
  i <- 1
  n <- length(cids)
  for (compound in cids) {
    tryCatch(
      {
        message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
        flush.console()

        url <- sprintf(url.base, compound)
        res <- fromJSON(url)
        cid <- c(cid, res[[1]][[1]]$CID)
        name <- c(name, unlist(res[[1]][[1]]$Synonym)[1])
        synonyms <- c(synonyms,
                      paste0(unlist(res[[1]][[1]]$Synonym), collapse = "; "))
        i <- i + 1
      },
      error = function(e){
        print(e)
      },
      finally = Sys.sleep(0.2)
    )
  }
  df <- data.frame(cid = cid, name = name, synonyms = synonyms,
                   stringsAsFactors = FALSE)
  return(df)
}

#' Get properties of drugs
#'
#' \code{GetPubchemPro} function retrieves the properties (InChIKey, Canonical
#' SMILES, and molecula formula) of drugs from PubChem database via
#' \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567}{PUG REST}.
#'
#' @param cids A vector of integer or character indicates the CIDs of drugs.
#'
#' @return A data frame contains 4 columns:
#' \itemize{
#'   \item \strong{CID} CID of drugs which is inputted to \code{cids} argument.
#'   \item \strong{InChIKey} Standard InChIKey of matched drugs.
#'   \item \strong{CanonicalSMILES} Standard Canonical SMILES of matched drugs.
#'   \item \strong{MolecularFormula} Molecular formula for matched drugs.
#' }
#' @export
#'
#' @examples
#' property <- GetPubchemPro(c(1,2,3,4))
GetPubchemPro <- function(cids){
  tryCatch({
  res <- NULL
  batch <- split(cids, ceiling(seq_along(cids)/100))

  for (i in 1:length(batch)) {
  temp <- NULL
  compound <- paste0(batch[[i]], collapse = ",")
  property <- paste0(c("InChIKey", "CanonicalSMILES", "MolecularFormula"),
                     collapse = ",")
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     compound, "/property/", property, "/CSV")
  temp <- read.csv(url, stringsAsFactors = FALSE)
  res <- rbind.data.frame(res, temp)
  }

  },
  error = function(e){
    print(e)
  },
  finally = Sys.sleep(0.2)
  )
  return(res)
}

GetPubPhase <- function(cids, quiet = TRUE){
  message("Getting clinical phases from PubChem.")
  # build container
  clinical_phase <- NULL
  # set indicator
  i <- 1
  n <- length(cids)

  for (compound in cids) {
    tryCatch({
      message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
      flush.console()

      url <- paste0('https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?',
                    'infmt=json&outfmt=json&query={"select":["cid","phase"],',
                    '"collection":"clinicaltrials",',
                    '"where":{"ands":[{"cid":"', compound, '"}]},',
                    '"order":["phase,desc"],',
                    '"start":1,"limit":10000}')
      res <- fromJSON(url)
      status <- res$SDQOutputSet[[1]][[1]]
      temp <- matrix(c(compound, 0), nrow = 1)
      colnames(temp) <- c("cid", "phase")

      if (status == 0) {
        # extract max clinical phase
        phase <- res$SDQOutputSet[[5]][[1]]
        if (length(phase) != 0) {
          temp <- res$SDQOutputSet[[5]][[1]][1,]
        }
      } else if (status != 0) {
        error <- res$SDQOutputSet[[1]][[2]]
        if (!quiet) {
          print( strsplit(error, ". ", fixed = T)[[1]][1] )
        }
      }

      clinical_phase <- rbind(clinical_phase, temp)
      i <- i + 1
    },
    error = function(e){
      if (!quiet) print(e)
    },
    finally = Sys.sleep(0.2)
    )
  }

  # clean
  gc()
  return(clinical_phase)
}

GetPubIDs <- fuction(cids) {
  url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                     "%s", "/xrefs/SourceName,RegistryID/JSON")
  i <- 1
  n <- length(cids)

  for (compound in cids) {
    message(round(i/n, 2)*100, "%", "completed", "\r", appendLF = FALSE)
    flush.console()

    tryCatch({
      url <- sprintf(url.base, compound)
      res <- fromJSON(url, flatten = TRUE)
      cid <- res[[1]][[1]]$CID
      ids <- res[[1]][[1]]$RegistryID[[1]]
      resource <- res[[1]][[1]]$SourceName[[1]]
    })

  }

}