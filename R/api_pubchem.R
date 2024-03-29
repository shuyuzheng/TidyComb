################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for retrieving or updating drug information from PubChem database.

#' Match CID according to other identifiers
#'
#' \code{GetCid} matches CID of drugs according to user-provided identifiers.
#'
#' \code{GetCid} queries PubChem database via
#' \href{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}{PUG REST}
#' to search matched CIDs of drugs according to other identifiers. Available
#' identifiers are: name, SMILES, InChIKey.
#'
#' Following is the \strong{"API USAGE POLICY"}:
#' Please note that PUG REST is not designed for very large volumes (millions)
#' of requests. We ask that any script or application not make more than 5
#' requests per second, in order to avoid overloading the PubChem servers. If
#' you have a large data set that you need to compute with, please contact us
#' for help on optimizing your task, as there are likely more efficient ways to
#' approach such bulk queries.
#'
#' @param ids A vector of characters contains the identifiers of drugs for
#' searching.
#'
#' @param type A character indicates the type of identifiers passed to \code{id}
#'   argument. Available types are:
#'   \itemize{
#'   \item \strong{smiles} Identifiers of drugs in the simplified molecular-input
#'   line-entry system (SMILES).
#'   \item \strong{inchikey} Standard InChIKey of the drugs.
#'   \item \strong{name} Name for drugs. This type could be used for searching
#'   synonyms, NCGC IDs, Chembl IDs, CAS or any other kind of identifiers.
#'   }
#'
#' @param quiet A logical value. If it is \code{TRUE}, the error message
#' during retrieving data will not show in console.
#'
#' @return A data frame contains two columns:
#' \itemize{
#'   \item \code{input_id} The identifiers input by user.
#'   \item \code{cid} The matched CID.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' GetCid("aspirin", "name", quiet = TRUE)
GetCid <- function(ids, type , quiet = TRUE) {
  message("Getting CIDs from PubChem...")
  types <- c("name", "smiles", "inchikey")

  if (!(type %in% types)) {
    stop("Invalid idtype specified, valiable idtypes are: ",
         paste(types, collapse = ", "))
  }

  cid <- integer()
  input_id <- character()
  url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/",
                       type, "/%s", "/cids/json")
  stepi <- 1
  n <- length(ids)
  for (id in ids) {

    temp <- NA
    tryCatch({

      url <- sprintf(url.base, utils::URLencode(id))

      doc <- jsonlite::fromJSON(url)
      rootNode <- names(doc)
      if (rootNode == "IdentifierList") {
        temp <- doc$`IdentifierList`$`CID`
      } else if (rootNode == "Fault") {
        fault <- doc$Fault$Details
        if (!quiet) {
          print( paste(id, fault[[1]], sep = ": ") )
        }
        temp <- NA
      }
    }, error = function(e) {
      if (!quiet) {
        print(e)
      }
    }, finally = Sys.sleep(0.2) # See usage policy.
    )
    cid <- c(cid, temp)
    input_id <- c(input_id, rep_len(id, length(temp)))

    message(round(stepi/n * 100), "%", "\r", appendLF = FALSE)
    utils::flush.console()

    stepi <- stepi + 1
  }
  df <- data.frame(input_id = input_id, cid = cid, stringsAsFactors = FALSE)
  return(df)
}

#' Get drug synonyms from PubChem
#'
#' \code{GetPubNames} queries PubChem database with drug CIDs via
#' (\href{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}{PUG REST} to
#' searching for synonyms of drug.
#'
#' @param cids A vector of characters or integers containing CID of drugs for
#' searching.
#' @param quiet A logical value. If it is \code{TRUE} the message for function
#' progress will be shutted down.
#'
#' @return A Data frame contains:
#' \itemize{
#'   \item \strong{cid} Inputted CIDs.
#'   \item \strong{name} The first name whose length is less than 40 characters
#'   in the synonyms list retrieved from PubChem.
#'   \item \strong{synonyms} Synonyms retrieved from PubChem
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' GetPubNames("2244")
GetPubNames <- function(cids, quiet = FALSE){
  if (!quiet) {
    message("Getting names from PubChem...")
  }
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
    tryCatch({
      if (!quiet) {
        message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
        utils::flush.console()
      }
      url <- sprintf(url.base, compound)
      res <- jsonlite::fromJSON(url)
      cid <- c(cid, res[[1]][[1]]$CID)
      names <- unlist(res[[1]][[1]]$Synonym)
      j <- 1
      dname <- paste0(rep(1, 41), collapse = "")
      while(nchar(dname) > 40 & j < length(names)){
        dname <- names[j]
        j <- j + 1
      }
      if (nchar(dname) > 40){
        dname <- names[1]
      }
      name <- c(name, dname)
      synonyms <- c(synonyms,
                    paste0(names, collapse = "; "))
    }, error = function(e){
      if (!quiet){
        print(e)
      }
      url <- sub("/synonyms", "", url)
      res <- jsonlite::fromJSON(url)
      cid <<- c(cid, res$PC_Compounds$id$id$cid)
      tmp <- res$PC_Compounds$props[[1]]
      names <- unique(tmp$value$sval[which(tmp$urn$label == "IUPAC Name")])
      j <- 0
      dname <- paste0(rep(1, 41), collapse = "")
      while(nchar(dname) > 40 & j < length(names)){
        j <- j + 1
        dname <- names[j]
      }
      if (nchar(dname) > 40){
        dname <- names[1]
      }
      name <<- c(name, dname)
      synonyms <<- c(synonyms,
                    paste0(names, collapse = "; "))
    }, finally = Sys.sleep(0.2)
    )
    i <- i + 1
  }
  df <- data.frame(cid = cid, name = name, synonyms = synonyms,
                   stringsAsFactors = FALSE)
  return(df)
}

#' Get properties of drugs
#'
#' \code{GetPubchemPro} function retrieves the properties (InChIKey, Canonical
#' SMILES, and molecula formula) of drugs from PubChem database via
#' \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567}{PUG REST},
#' accordint to CIDs.
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
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @examples
#' property <- GetPubchemPro(c(1,2,3,4))
GetPubchemPro <- function(cids) {
  message("Getting drug properties from PubChem...")

  res <- NULL
  batch <- split(cids, ceiling(seq_along(cids)/100))

  for (i in 1:length(batch)) {
    tryCatch({
      temp <- NULL
      compound <- paste0(batch[[i]], collapse = ",")
      property <- paste0(c("InChIKey", "CanonicalSMILES", "MolecularFormula"),
                         collapse = ",")
      url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                         compound, "/property/", property, "/CSV")
      temp <- utils::read.csv(url, stringsAsFactors = FALSE)
      res <- rbind.data.frame(res, temp)
    }, error = function(e){
        print(e)
    }, finally = Sys.sleep(0.2)
    )
  }
      colnames(res) <- c("cid", "inchikey", "smiles", "molecular_formula")
      return(res)
}

#' Get max clinical trial phase from Pubchem
#'
#' \code{GetPubPhase} function retrieves the max clinical trial phase of drugs
#' from PubChem database via
#' \href{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865567}{PUG REST},
#' accordint to CIDs.
#'
#' @param cids A vector of integer or character indicates the CIDs of drugs.
#'
#' @param quiet A logical value. If it is \code{TRUE}, the error message will
#' shut down during retrieving process.
#'
#' @return A data frame contains 4 columns:
#' \itemize{
#'   \item \strong{CID} CID of drugs which is inputted to \code{cids} argument.
#'   \item \strong{phase} Max clinical trial phase of matched drugs.
#' }
#' @export
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @examples
#' clinical.phase <- GetPubPhase(c(1,2,3,4, 2244))
GetPubPhase <- function(cids, quiet = TRUE) {
  message("Getting clinical phases from PubChem...")
  # build container
  clinical_phase <- data.frame(cid = integer(),
                               phase = integer(),
                               stringsAsFactors = FALSE)
  # set indicator
  i <- 1
  n <- length(cids)

  for (compound in cids) {

    message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
    utils::flush.console()

    tryCatch({
      url <- paste0('https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?',
                    'infmt=json&outfmt=json&query={"select":["cid","phase"],',
                    '"collection":"clinicaltrials",',
                    '"where":{"ands":[{"cid":"', compound, '"}]},',
                    '"order":["phase,desc"],',
                    '"start":1,"limit":10000}')
      res <- jsonlite::fromJSON(url)
      status <- res$SDQOutputSet$status$code

      if (status == 0) {
        # extract max clinical phase
        phase <- unique(res$SDQOutputSet$rows[[1]]$phase)
        phase <- sub("N/A", "", phase)
        phase <- sub("^.*/", "", phase)
        phase <- sub("Early Phase 1", "0", phase)
        phase <- sub("Phase ", "", phase)
        phase <- unique(phase)
        if (length(phase) == 0) {
          phase <- ""
        } else if (length(phase) == 1) {
          phase <- phase
        } else if (length(phase) > 1) {
          phase <- max(as.integer(phase), na.rm = TRUE)
        }
      } else {
        error <- res$SDQOutputSet[[1]][[2]]
        phase <- ""
        if (!quiet) {
          print( strsplit(error, ". ", fixed = T)[[1]][1])
        }
      }
    }, error = function(e){
      phase <<- ""
      if (!quiet) {
        print(e)
      }
    }, finally = {
      Sys.sleep(0.2)
      temp <- data.frame(cid = compound,
                         phase = phase,
                         stringsAsFactors = FALSE)
    }
    )

    clinical_phase <- rbind.data.frame(clinical_phase, temp)
    i <- i + 1
    temp <- NULL
  }

  return(clinical_phase)
}

#' Get drug synonyms from PubChem by Drug names
#'
#' \code{GetPubSynonymFromName} queries PubChem database with drug names via
#' (\href{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}{PUG REST} and
#' retrieves synonyms of drugs.
#'
#' @param names A vector of characters containing names of drugs for searching.
#'
#' @return A Data frame contains:
#' \itemize{
#'   \item \strong{input_name} The name inputted by user.
#'   \item \strong{synonyms} Synonyms retrieved from PubChem. Different synonyms
#'   are separated by ";".
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' names <- c("Aspirin", "5-FU")
#' synonyms <- GetPubSynonymFromName(names)
#'

GetPubSynonymFromName <- function(names) {
  message("Getting names from PubChem...")
  url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                     "%s", "/synonyms/JSON")

  # Build containers
  input_name <- NULL
  synonyms <- NULL
  n <- length(names)
  for (i in 1:n) {
    message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
    utils::flush.console()
    tryCatch({
      url <- sprintf(url.base, names[i])
      res <- jsonlite::fromJSON(url)
      input_name <- c(input_name, names[i])
      synonyms <- c(synonyms,
                    paste0(unlist(res[[1]][[1]]$Synonym), collapse = "; "))
    }, error = function(e){
      print(e)
    }, finally = Sys.sleep(0.2)
    )
  }
  df <- data.frame(input_name = input_name, synonyms = synonyms,
                   stringsAsFactors = FALSE)
  return(df)
}

#' Get "Preferred Compound" of "Non-live" compounds in PubChem
#'
#' Fuction \code{UpdateCid} takes CIDs as input. If inputted compounds was
#' labeled as "Non-live" in PubChem, the CID of the "preferred Compound" will
#' be returned. If the inputted CID is still "live" in PubChem, \code{NA} will
#' be returned.
#'
#' @param cids a vector contains the CIDs of drugs on which test will be done.
#'
#' @return a data frame contains 2 columns:
#'   \itemize{
#'     \item \strong{old_cid} The inputted CIDs.
#'     \item \strong{new_cid} If the tested compound is "Non-live", it is the
#'     CID of "preferred compound". If the tested compound is still "live", it
#'     is \code{NA}.
#'   }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
#'
#' @examples
#' res <- UpdateCid(1:10)
UpdateCid <- function(cids) {
  url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/",
                     "compound/%s/JSON")
  output <- NULL
  n <- length(cids)

    for (i in 1:n) {
      message(round(i/n, 2)*100, "% completed", "\r", appendLF = FALSE)
      utils::flush.console()
      tmp <- data.frame(old_cid = cids[i], new_cid = NA, stringsAsFactors = FALSE)
      new_cid <- NA
      tryCatch({
        url <- sprintf(url.base, cids[i])
        res <- jsonlite::fromJSON(url)
        if ("Preferred Compound" %in% res$Record$Section$TOCHeading) {
          new_cid <- res$Record$Section$Information[[1]]$Value$Number
        }
      }, error = function(e) {
        print(paste0("Error on ", cids[i], ": ",  e))
        new_cid <- NA
      }, finally = Sys.sleep(0.2)
      )

      tmp$new_cid <- new_cid
      output <- rbind.data.frame(output, tmp)
    }

  return(output)
}
# GetPubIDs <- fuction(cids) {
#   url.base <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
#                      "%s", "/xrefs/SourceName,RegistryID/JSON")
#   i <- 1
#   n <- length(cids)
#
#   for (compound in cids) {
#     message(round(i/n, 2)*100, "%", "completed", "\r", appendLF = FALSE)
#     utils::flush.console()
#
#     tryCatch({
#       url <- sprintf(url.base, compound)
#       res <- fromJSON(url, flatten = TRUE)
#       cid <- res[[1]][[1]]$CID
#       ids <- res[[1]][[1]]$RegistryID[[1]]
#       resource <- res[[1]][[1]]$SourceName[[1]]
#     })
#
#   }
#
# }

# GetPugKEGGXML <- function(cids) {
#   query <- XML::xmlParse(system.file('api',
#                                     "pugquery.xml",
#                                     package='TidyComb'))
#   parent <- XML::getNodeSet(query, "//PCT-ID-List_uids")[[1]]
#   n <- lapply(cids, function(x) {XML::newXMLNode("PCT-ID-List_uids_E", x)})
#   addChildren(parent, kids = n)
#   query <- saveXML(query)
#   query <- gsub("\\n", "", query)
#   query <- gsub('<\\?xml version="1\\.0"\\?>', '', query)
#   return(query)
# }
#
# GetPollBody <- function(reqid) {
#   doc <- XML::xmlParse(system.file('api',
#                                    "pugpoll.xml",
#                                    package='TidyComb'))
#   node <- XML::getNodeSet(doc, "//PCT-Request_reqid")[[1]]
#   XML::xmlValue(node) <- reqid
#   doc <- saveXML(doc)
#   doc <- gsub('<\\?xml version="1\\.0"\\?>', '', doc)
#   doc <- gsub("\\n", "", doc)
#   return(doc)
# }


# h <- RCurl::basicTextGatherer()
# RCurl::curlPerform(url = "http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi",
#             postfields = query,
#             writefunction = h$update)
# xml <- xmlTreeParse(h$value(), asText = TRUE, asTree = TRUE)
# root <- xmlRoot(xml)
#
#
#
#
# PollPug <- function(reqid) {
#   root <- NA
#   pstring <- GetPollXml(reqid)
#   reqid <- NA
#   while(TRUE) {
#     h = RCurl::basicTextGatherer()
#     RCurl::curlPerform(url = 'http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi',
#                       postfields = pstring,
#                       writefunction = h$update)
#     ## see if we got a waiting response
#     root <- xmlRoot(xmlTreeParse(h$value(), asText=TRUE, asTree=TRUE))
#     reqid <- xmlElementsByTagName(root, 'PCT-Waiting', recursive=TRUE)
#     if (length(reqid) != 0) next
#     break
#   }
#   return(root)
# }
#
# .get.aid.by.cid.old <- function(cid, type='raw', quiet=TRUE) {
#
#   if (!(type %in% c('tested','active','inactive','discrepant','raw')))
#     stop("Invalid type specified")
#
#   url <- "http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
#
#   ## perform query
#   qstring <- gsub("\\n", "", sprintf(.queryString, cid))
#   h = basicTextGatherer()
#   curlPerform(url = 'http://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi',
#               postfields = qstring,
#               writefunction = h$update)
#
#   ## extract query id
#   xml <- xmlTreeParse(h$value(), asText=TRUE, asTree=TRUE)
#   root <- xmlRoot(xml)
#   reqid <- xmlElementsByTagName(root, 'PCT-Waiting_reqid', recursive=TRUE)
#   if (length(reqid) != 1) {
#     if (!quiet) warning("Malformed request id document")
#     return(NULL)
#   }
#   reqid <- xmlValue(reqid[[1]])
#
#   ## start polling
#   if (!quiet) cat("Starting polling using reqid:", reqid, "\n")
#   root <- .poll.pubchem(reqid)
#
#   ## OK, got the link to our result
#   link <- xmlElementsByTagName(root, 'PCT-Download-URL_url', recursive=TRUE)
#   if (length(link) != 1) {
#     if (!quiet) warning("Polling finished but no download URL")
#     return(NULL)
#   }
#   link <- xmlValue(link[[1]])
#   if (!quiet) cat("Got link to download:", link, "\n")
#
#   ## OK, get data file
#   tmpdest <- tempfile(pattern = 'abyc')
#   tmpdest <- paste(tmpdest, '.gz', sep='', collapse='')
#   status <- try(download.file(link,
#                               destfile=tmpdest,
#                               method='internal',
#                               mode='wb', quiet=TRUE),
#                 silent=TRUE)
#   if (class(status) == 'try-error') {
#     if (!quiet) warning(status)
#     return(NULL)
#   }
#
#   ## OK, load the data
#   dat <- utils::read.csv(tmpdest,header=TRUE,fill=TRUE,row.names=NULL)
#   unlink(tmpdest)
#
#   valid.rows <- grep("^[[:digit:]]*$", dat[,1])
#   dat <- dat[valid.rows,c(1,3,4,5)]
#   row.names(dat) <- 1:nrow(dat)
#   names(dat) <- c('aid', 'active', 'inactive', 'tested')
#   ret <- dat
#
#   type <- type[1]
#   switch(type,
#          active = dat[dat$active == 1,1],
#          inactive = dat[dat$inactive == 1,1],
#          tested = dat[,1],
#          raw = ret[,-5])
# }