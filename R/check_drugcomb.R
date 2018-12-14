#' Check current status of DrugComb
#'
#' Get the updated csv files from DrugComb database:
#' cell_line: cellosaurus_accession, id
#' select * from drug into outfile '/var/lib/mysql-files/dr.csv' fields terminated by ',' enclosed by '"' lines terminated by '\n';
#' drug: id, cid
#' tissue: id, name (*)
#' disease: name, id (*)
#' study: id
#'
CheckCell <- function(test) {
  message("Checking Cell lines...")
  exist <- read.csv(system.file("extdata", "cell_line.csv",
                                package = "TidyComb"),
                    header = FALSE,
                    col.names = c("cellosaurus_accession", "id"),
                    stringsAsFactors = FALSE)
  n <- nrow(exist)
  old <- exist[match(test, exist$cellosaurus_accession), ] %>% na.omit()
  new <- test[!test %in% exist$cellosaurus_accession]
  message("DrugComb has archived ", n, "cell lines.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " of checked cell line(s) have/has been in DrugComb: ",
          paste(old$cellosaurus_accession, collapse = ", "), "\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " of checked cell line(s) are/is not in DrugComb:",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}

CheckTissue <- function(test) {
  message("Checking tissues...")
  exist <- read.csv(system.file("extdata", "tissue.csv",
                                package = "TidyComb"),
                    header = FALSE,
                    col.names = c("id", "name"),
                    stringsAsFactors = FALSE)
  n <- nrow(exist)
  old <- exist[match(test, exist$name), ] %>% na.omit()
  new <- test[!test %in% exist$name]
  message("DrugComb has archived ", n, "tissues.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " of checked tissue(s) have/has been in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " of checked tissue(s) are/is not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}

CheckDisease <- function(test) {
  message("Checking  diseases...")
  exist <- read.csv(system.file("extdata", "disease.csv",
                                package = "TidyComb"),
                    header = FALSE,
                    col.names = c("name", "id"),
                    stringsAsFactors = FALSE)
  n <- nrow(exist)
  old <- exist[match(test, exist$id), ] %>% na.omit()
  new <- test[!test %in% exist$id]
  message("DrugComb has archived ", n, "diseases.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " of checked disease(s) have/has been in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " of checked disease(s) are/is not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}

CheckDrug <- function(cids) {
  message("Checking drugs...")
  exist <- read.csv(system.file("extdata", "drug.csv",
                                package = "TidyComb"),
                    header = FALSE,
                    col.names = c("id", "cid"),
                    stringsAsFactors = FALSE)
  n <- nrow(exist)
  old <- exist[match(cids, exist$cid), ] %>% na.omit()
  new <- cids[!cids %in% exist$cid]
  message("DrugComb has archived ", n, "drugs.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " of checked drug(s) have/has been in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " of checked drug(s) are/is not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}