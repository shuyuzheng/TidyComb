
library(XML)

# 1. Check data version ---------------------------------------------------

CellVersion <- function(){
  #build handler
  startElement = function(ctxt, name, attrs, ...) {
    if (name == "release") {
      print(attrs)
      xmlStopParser(ctxt)
    }
  }
  class(startElement) = "XMLParserContextFunction" 
  
  # check online c
  local <- file.path("raw", "cellosaurus.xml")
  url <- "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml"
  print("Online version is:")
  invisible(xmlEventParse(url, handlers = list(startElement = startElement)))
  print("Local version is:")
  invisible(xmlEventParse(local, handlers = list(startElement = startElement)))
}

# 2. Update local file if the online version updated. ---------------------

download.file("ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml",
              file.path("raw", "cellosaurus.xml"))

# 3. Get cell information, according to cell name -------------------------
#' Extract cell informations from cellosaurus's data
#' 
#' 
GetCell <- function(cells, file, quiet = TRUE) {
  
  # Input: character vector of cell lines' name
  # Output: data.frame with matched cellosaurus accession ID, tissue and disease

  tryCatch({
    
  # parsing file
  file <- file.path("raw", "cellosaurus.xml")
  doc <- xmlInternalTreeParse(file = file)
  
  # Building container
  info <- NULL
  
  stepi <- 0
  n <- length(cells)

  for (cell in cells) {
    # set a pregressing index
    message(round(stepi/n * 100), "%", " is completed",
            "\r", appendLF = FALSE)
    flush.console()
    

    # building containers
    synonyms <- character()
    disease <- character()
    disease_id <- character()
    accession <- character()
    name <- character()  
    
    # catch node with the "cell name"
    node <- xpathApply(doc, path = paste0("//name[text()=\'", cell,
                                          "\']/../.."),xmlToList)
    
    node.info <- NULL
    for (i in length(node)) {
      new_item <- NA
      # ger synonyms
      new_item <- unlist(sapply(node[[i]]$`name-list`, "[", "text"))
        if (is.null(new_item)) {
            new_item <- NA
          }else{
            new_item <- paste(new_item[-1], collapse = ";")
          }
        
      synonyms <- c(synonyms, new_item)
      new_item <- NA
      
      # get accession ID
      new_item <- unlist(sapply(node[[i]]$"accession-list","[", "text"))
        if (is.null(new_item)) {
         new_item <- NA
        } 
      accession <- c(accession, new_item)
      new_item <- NA

      # get disease type
      new_item <- unlist(node[[i]]$"disease-list")
      if (length(new_item) == 0) {
        new_item <- NA
      }
      disease <- c(disease, new_item["cv-term.text"])
      disease_id <- c(disease_id, new_item["cv-term..attrs.accession"])
      new_item <- NA
      
      # record the cell name
      name <- c(name, cell)
      
      # combine the data of one node
      node.info <- data.frame(name, synonyms, cellosaurus_accession = accession,
                         disease, disease_id, stringsAsFactors = FALSE)
    }
    
    # combine the data of one cell
    info <- rbind.data.frame(info, node.info)
    
    stepi <- stepi + 1
    }
  },
  error = function(e) {
    if (!quiet) {
      print(e)
    }
  
  })

  # clean
  
  gc()
  return(info)
}
