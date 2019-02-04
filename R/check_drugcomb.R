# TidyComb
# Functions for checking current status of DrugComb
# Copyright: Shuyu Zheng
#
# All functions in this section are based on an up-to-date csv files from
# DrugComb database. It consists of :
#
# cell_line: cellosaurus_accession, id
# drug: id, cid
# tissue: id, name (*)
# disease: name, id (*)
# study: id
#
# SQL command for querying data from DrugComb database:
#
# select [variables] from [table] into outfile '/var/lib/mysql-files/dr.csv'
# fields terminated by ',' enclosed by '"' lines terminated by '\n';

#' Check cells archive in DrugComb
#'
#' Function \code{CheckCell} compares cellosaurus accessions passed to
#' \code{test} argument with cell lines archived in DrugComb (extracted to
#' \code{\link{drugcomb}} R data). It will return the total number of cell lines
#' have been archived in DrugComb as well as which cell lines are already or not
#' yet archived in DrugComb.
#'
#' @param test A character vector. It contains \emph{cellosaurus accessions}
#' of cell lines which would be checked.
#'
#' @return A list contains 3 elements:
#' \itemize{
#'   \item \strong{n} It indicates that there are totally n cell lines archived
#'   in DrugComb database.
#'   \item \strong{old} It lists out the \emph{cellosaurus accessions} of cell
#'   lines that \emph{have been archived} in DrugComb.
#'   \item \strong{new} It lists out the \emph{cellosaurus accessions} of cell
#'   lines that \emph{are not yet} in DrugComb.
#' }
#'
#' @export
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

#' Check tissues archive in DrugComb
#'
#' Function \code{CheckTissue} compares tissue names passed to \code{test}
#' argument with tissues archived in DrugComb (extracted to
#' \code{\link{drugcomb}} R data). It will return the total number of tissues
#' have been archived in DrugComb as well as which tissues are already or not
#' yet archived in DrugComb.
#'
#' @param test A character vector. It contains names of tissues which would be
#' checked.
#'
#' @return A list contains 3 elements:
#' \itemize{
#'   \item \strong{n} It indicates that there are totally n tissues archived
#'   in DrugComb database.
#'   \item \strong{old} It lists out the names of tissues that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new} It lists out the names of tissues that
#'   \emph{are not yet} in DrugComb.
#' }
#'
#' @export
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

#' Check diseases archive in DrugComb
#'
#' Function \code{CheckDisease} compares disease IDs passed to \code{test}
#' argument with diseases archived in DrugComb (extracted to
#' \code{\link{drugcomb}} R data). It will return the total number of diseases
#' have been archived in DrugComb as well as which diseases are already or not
#' yet archived in DrugComb.
#'
#' @param test A character vector. It contains IDs of diseases which would
#' be checked.
#'
#' @return A list contains 3 elements:
#' \itemize{
#'   \item \strong{n} It indicates that there are totally n diseases archived
#'   in DrugComb database.
#'   \item \strong{old} It lists out the ID of diseases that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new} It lists out the ID of diseases that \emph{are not yet}
#'   in DrugComb.
#' }
#'
#' @export
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

#' Check drug archives in DrugComb
#'
#' Function \code{CheckDrug} compares drug CIDs passed to \code{cids} argument
#' with CIDs archived in DrugComb (extracted to \code{\link{drugcomb}} R data).
#' It will return the total number of drugs have been archived in DrugComb as
#' well as which drugs are already or not yet archived in DrugComb.
#'
#' @param cids A character vector. It contains CIDs of drugs which would be
#' checked.
#'
#' @return A list contains 3 elements:
#' \itemize{
#'   \item \strong{n} It indicates that there are totally n drugs archived in
#'   DrugComb database.
#'   \item \strong{old} It lists out the CID of drugs that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new} It lists out the CID of drugs that \emph{are not yet}
#'   in DrugComb.
#' }
#'
#' @export
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