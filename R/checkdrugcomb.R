#' Check current status of DrugComb
#'
#' Get the updated csv files from DrugComb database:
#' cell_line: synonyms, cellosaurus_accession, id
#' select * from drug into outfile '/var/lib/mysql-files/dr.csv' fields terminated by ',' enclosed by '"' lines terminated by '\n';
#' drug: id, cid
#' tissue: id, name (*)
#' disease: name, id (*)
#' study: id
#'
CheckCell <- function(test) {
  exist <- read.csv(system.file("extdata", "cell_line.csv",
                                package = "TidyComb"),
                    col.names = c("cellosaurus_accession", "id"))
  n <- nrow(exist)
  old <- exist[match(test, exist$cellosaurus_accession), ] %>% na.omit()
  new <- test[!test %in% exist$cellosaurus_accession]
  message("There are ", n, " cell lines in DrugComb.", "\n",
      "Cell lines are in DrugComb: ",
      paste(old$cellosaurus_accession, collapse = ", "), "\n",
      "Cell lines are not in DrugComb: ",
      paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}

CheckTissue <- function(test) {
  exist <- read.csv(system.file("exdata", "tissue.csv",
                                package = "TidyComb"),
                    col.names = c("id", "name"))
  n <- nrow(exist)
  old <- exist[match(test, exist$name), ] %>% na.omit()
  new <- test[!test %in% exist$name]
  message("There are ", n, " tissues in DrugComb.", "\n",
          "Tissues are in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          "Tissues are not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(old = old, new = new))
}

CheckDisease <- function(test) {
  exist <- read.csv(system.file("exdata", "disease.csv",
                                package = "TidyComb"),
                    col.names = c("name", "id"))
  n <- nrow(exist)
  old <- exist[match(test, exist$id), ] %>% na.omit()
  new <- test[!test %in% exist$id]
  message("There are ", n, " disaeses in DrugComb.", "\n",
          "Diseases are in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          "Diseases are not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}

CheckDrug <- function(test) {
  exist <- read.csv(system.file("exdata", "drug.csv",
                                package = "TidyComb"),
                    col.names = c("id", "cid"))
  n <- nrow(exist)
  old <- exist[match(test, exist$cid), ] %>% na.omit()
  new <- test[!test %in% exist$cid]
  message("There are ", n, " drugs in DrugComb.", "\n",
          "Drugs are in DrugComb: ",
          paste(old$name, collapse = ", "), "\n",
          "Drugs are not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(old = old, new = new))
}