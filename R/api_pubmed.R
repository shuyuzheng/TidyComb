get.pubmed <- function(pmids, tool = "R", email = NULL ) {
  # Input:
  #    1. pmid: PubMed ID of publications
  #    2. email: user's e-mail address
  #    3. tool: the name of application, default is "R"
  # Output: data.frame with matched publication date, first author, title,
  #    journal, and doi of the publication.
  #
  # API Documentation: https://www.ncbi.nlm.nih.gov/pmc/tools/developers/
  #
  # Note about API usage:
  # The use of our APIs is entirely free, and doesn't require an API key, but
  # we ask that you please:
  # * Do not make concurrent requests, even at off-peak times;
  # * Include two parameters that help to identify your service or application
  # to our servers:
  #    1. tool: should be the name of the application, as a string value with
  # no internal spaces.
  #    2. email: should be the e-mail address of the maintainer of the tool,
  # and should be a valid e-mail address.
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

