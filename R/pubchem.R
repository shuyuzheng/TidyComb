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


# ----------2. get.info function get drug information according to CIDs.

GetDrug <- function(cids, quiet = TRUE){

  # Input: character vector of compound CIDs
  # Output: data.frame with matched names, InChI Keys, SMILES, molecula
  #         formulas.
  #
  # API Documentation: https://pubchemdocs.ncbi.nlm.nih.gov/pug-view

  tryCatch({
  # build containers
  inchikey <- character()
  smiles <- character()
  molecula_formula <- character()
  clinical_phase <- character()
  name <- character()
  cid <- integer()
  for (compound in cids) {
    message(compound, "\r", appendLF = FALSE)
    flush.console()
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data",
                    "/compound/", compound, "/XML")
      curlHandle <- getCurlHandle()
      res <- dynCurlReader()
      curlPerform(url = url, curl = curlHandle, writefunction = res$update)
      doc <- htmlTreeParse(res$value(), useInternalNodes = TRUE)
      body <- xpathApply(doc, "//body", xmlChildren)
      new_item <- NA
      # extract informations
    if (names(body[[1]]) == "record") {
      # extract inchikey
      new_item <- xpathApply(doc,
                    "//tocheading[contains(text(), 'InChI Key')]",
                    function(x){
                      xpathSApply(x, "..//stringvalue", xmlValue)
                    })
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
      inchikey <- c(inchikey, new_item)

      # extract smiles
      new_item <- xpathApply(doc,
                    "//tocheading[contains(text(), 'Canonical SMILES')]",
                    function(x){
                      xpathSApply(x, "..//stringvalue", xmlValue)
                    })
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
     smiles <- c(smiles, new_item)

      # extract molecula formula
      new_item <- xpathApply(doc,
                    "//tocheading[contains(text(), 'Molecular Formula')]",
                    function(x) {
                      xpathSApply(x, "..//stringvalue", xmlValue)
                    })
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
      molecula_formula <- c(molecula_formula, new_item)

      # extract names
      new_item <- xpathApply(doc,
                  "//description[contains(text(), 'Primary Identifying Name')]",
                  function(x){
                    xpathApply(x, "..//information/stringvalue", xmlValue)
                  })
      new_item <- unlist(new_item)
      if (is.null(new_item)) {
        new_item <- NA
      }
      name <- c(name, unlist(new_item))
    } else if (names(body[[1]]) == "fault") {
      fault <- xpathApply(doc, "//message", xmlValue)
      if (!quiet) {
        print( paste(compound, fault[[1]], sep = ": ") )
      }
      inchikey <- c(inchikey, NA)
      smiles <- c(smiles, NA)
      molecula_formula <- c(molecula_formula, NA)
      name <- c(name, NA)
      clinical_phase <- c(clinical_phase, NA)
     }
     cid <- c(cid, compound)
    }
  info <- cbind(name, inchikey, smiles, cid, molecula_formula, clinical_phase)
  },
  error = function(e){
    if (!quiet) {
      print(e)
    }
  },
  finally = Sys.sleep(0.2)
  )

  # clean
  remove(curlHandle)
  gc()
  return(info)
}

GetSid <- function(cids) {
  url.base <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/sids/JSON"
  df <- data.frame()
  for (cid in cids) {
    message(cid, "\r", appendLF = FALSE)
    flush.console()

    url <- sprintf(url.base, cid)
    temp <- as.data.frame(fromJSON(url))
    df <- rbind.data.frame(df, temp)
    temp <- NA
  }
  return(df)
}
GetInchikey <- function(cids, quiet = TRUE) {
  tryCatch({
    # build containers
    url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                      "%s/property/InChIKey/JSON")
    df <- data.frame()
    for (cid in cids) {
      message(cid, "\r", appendLF = FALSE)
      flush.console()

      url <- sprintf(url.base, cid)
      temp <- as.data.frame(fromJSON(url))
      df <- rbind.data.frame(df, temp)
      temp <- NA
    }
  },
  error = function(e) {
    print(e)
    temp <- data.frame(cid = cid, inchikey = NA)
    df <- rbind.data.frame(df, temp)
    temp <- NA
  },
  finally = Sys.sleep(0.2)
  )

  # clean
  gc()
  return(df)
}

GetClinicalPhase <- function(cids, quiet = TRUE){

  # Input: character vector of compound CIDs
  # Output: data.frame with matched max clinical phase.
  #
  # API: "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi"
  #
  # USAGE POLICY: Please note that PUG View is not designed for very large
  # volumes (millions) of requests. We ask that any script or application not
  # make more than 5 requests per second, in order to avoid overloading the
  # PubChem servers. To check additional request volume limitations, please read
  # this document. If you have a large data set that you need to compute with,
  # please contact us for help on optimizing your task, as there are likely more
  # efficient ways to approach such bulk queries.

  tryCatch({
  # build containers
  clinical_phase <- data.frame()
  for (compound in cids) {
    message(compound, "\r", appendLF = FALSE)
    flush.console()
      url <- paste0('https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?',
                    'infmt=json&outfmt=json&query={"select":["cid","phase"],',
                    '"collection":"clinicaltrials",',
                    '"where":{"ands":[{"cid":"', compound, '"}]},',
                    '"order":["phase,desc"],',
                    '"start":1,"limit":10000}')
      res <- fromJSON(url)
      status <- res$SDQOutputSet[[1]][[1]]

    if (status == 0) {
      # extract max clinical phase
      phase <- res$SDQOutputSet[[5]][[1]]
      if (length(phase) == 0) {
        new_item <- data.frame(phase = NA, cid = compound)
        clinical_phase <- rbind.data.frame(clinical_phase, new_item)
      } else{
        new_item <- res$SDQOutputSet[[5]][[1]][1,]
        clinical_phase <- rbind.data.frame(clinical_phase, new_item)
      }

   } else if (status != 0) {
      error <- res$SDQOutputSet[[1]][[2]]
      if (!quiet) {
        print( strsplit(error, ". ", fixed = T)[[1]][1] )
      }
      new_item <- data.frame(phase = NA, cid = compound)
      clinical_phase <- rbind.data.frame(clinical_phase, new_item)
     }
    }
  },
  error = function(e){
    print(e)
  },
  finally = Sys.sleep(0.2)
  )

  # clean
  gc()
  return(clinical_phase)
}

GetCidFromAid <- function(aids, quiet = TRUE){
   # Input: character vector of compound CIDs
  # Output: data.frame with matched names, InChI Keys, SMILES, molecula
  #         formulas.
  #
  # API Documentation: https://pubchemdocs.ncbi.nlm.nih.gov/rest/pug
  # API Documentation: https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
  #
  # USAGE POLICY: Please note that PUG REST is not designed for very large
  # volumes (millions) of requests. We ask that any script or application not
  # make more than 5 requests per second, in order to avoid overloading the
  # PubChem servers. If you have a large data set that you need to compute
  # with, please contact us for help on optimizing your task, as there are
  # likely more efficient ways to approach such bulk queries.

  tryCatch({
  # build containers
  res <- data.frame()
  for (assay in aids) {
    message(assay, "\r", appendLF = FALSE)
    flush.console()
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/",
                    assay,"/cids/json")
      doc <- fromJSON(url)
      new_item <- NA
      # extract informations
    if (names(doc) == "InformationList") {
      aid <- doc$InformationList[[1]][1]
      cid <- unlist(doc$InformationList[[1]][2])
      new_item <- cbind(aid, cid)
    } else if (names(doc) == "Fault") {
      fault <- doc$Fault[[2]]
      if (!quiet) {
        print( paste(assay, fault, sep = ": ") )
      }
     }
  }
  res <- rbind(res, new_item)
  },
  error = function(e){
    print(e)
  },
  finally = Sys.sleep(0.2)
  )

  # clean
  gc()
  return(res)
}
