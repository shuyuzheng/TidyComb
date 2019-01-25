# Check current status of DrugComb
#
# Get the updated csv files from DrugComb database:
# cell_line: cellosaurus_accession, id
# select * from drug into outfile '/var/lib/mysql-files/dr.csv' fields terminated by ',' enclosed by '"' lines terminated by '\n';
# drug: id, cid
# tissue: id, name (*)
# disease: name, id (*)
# study: id
#
CheckCell <- function(test) {
  message("Checking Cell lines...")
  exist <- drugcomb$cell_line
  n <- nrow(exist)
  old <- stats::na.omit(exist[match(test, exist$cellosaurus_accession), ])
  new <- test[!test %in% exist$cellosaurus_accession]
  message("DrugComb has archived ", n, " cell lines.\n",
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
  exist <- drugcomb$tissue
  n <- nrow(exist)
  old <- stats::na.omit(exist[match(test, exist$name), ])
  new <- test[!test %in% exist$name]
  message("DrugComb has archived ", n, " tissues.\n",
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
  exist <- drugcomb$disease
  n <- nrow(exist)
  old <- stats::na.omit(exist[match(test, exist$id), ])
  new <- test[!test %in% exist$id]
  message("DrugComb has archived ", n, " diseases.\n",
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
  exist <- drugcomb$drug
  n <- nrow(exist)
  old <- stats::na.omit(exist[match(cids, exist$cid), ])
  new <- cids[!cids %in% exist$cid]
  message("DrugComb has archived ", n, " drugs.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " of checked drug(s) have/has been in DrugComb: ",
          paste(old$cid, collapse = ", "), "\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " of checked drug(s) are/is not in DrugComb: ",
          paste0(new, collapse = ", "))
  return(list(n = n, old = old, new = new))
}