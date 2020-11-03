################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

# TidyComb
# Functions for checking current status of DrugComb
#
# Fuctions on this page:
# CheckCell: Check whether input cell lines are archived in DrugComb
# CheckTissue: Check whether input tissues are archived in DrugComb
# CheckDisease: Check whether input diseases are archived in DrugComb
# CheckDrug: Check whether input drugs are archived in DrugComb
#
# All functions in this section are based on up-to-date csv files from
# DrugComb database. It consists of :
#
# cell_line: cellosaurus_accession, id
# drug: id, cid
# tissue: id, name
# disease: name, id
# study: id
#
# SQL command for querying data from DrugComb database:
#
# select [variables] from [table] into outfile '/var/lib/mysql-files/dr.csv'
# fields terminated by ',' enclosed by '"' lines terminated by '\n';

#' Check whether input cell lines are archived in DrugComb
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
#'   \item \strong{n}: It indicates that there are totally n cell lines archived
#'   in DrugComb database.
#'   \item \strong{old}: It lists the \emph{cellosaurus accessions} of cell
#'   lines that \emph{have been archived} in DrugComb.
#'   \item \strong{new}: It lists the \emph{cellosaurus accessions} of cell
#'   lines that \emph{are not yet} in DrugComb.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
CheckCell <- function(test) {
  message("Checking Cell lines...")
  exist <- drugcomb$cell_line
  n <- nrow(exist)
  test <- stats::na.omit(test)
  old <- stats::na.omit(exist[match(test, exist$cellosaurus_accession), ])
  new <- test[!test %in% exist$cellosaurus_accession]
  message("DrugComb has archived ", n, " cell lines.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " inputted cell line(s) have/has been in DrugComb.\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " inputted cell line(s) are/is not in DrugComb.")
  return(list(n = n, old = old, new = new))
}

#' Check whether input tissues are archived in DrugComb
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
#'   \item \strong{n}: It indicates that there are totally n tissues archived
#'   in DrugComb database.
#'   \item \strong{old}: It lists the names and IDs of tissues that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new}: It lists the names of tissues that
#'   \emph{are not yet} in DrugComb.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
CheckTissue <- function(test) {
  message("Checking tissues...")
  if (sum(is.na(test)) > 0) {
    message("There are missing values in input data. Please check the ",
    "'cell_line' table and  manually fix missing tissues or remove missin",
    "values from input data with function 'na.omit'.")
    return(NULL)
  }
  if ("central_nervous_system" %in% test){
    message("Please check the cell lines from 'central_nevours_system' to specify
         which of the following subtypes they are from: 'brain',
         'nervours_system', or 'spinal_cord_or_other_CNS'.")
    return(NULL)
  }
  exist <- drugcomb$tissue
  n <- nrow(exist)
  test <- stats::na.omit(test)
  old <- stats::na.omit(exist[match(test, exist$tname), ])
  new <- test[!test %in% exist$tname]
  new <- data.frame(id = seq(n + 1, length.out = length(new)),
                    tname = new,
                    stringsAsFactors = FALSE)
  message("DrugComb has archived ", n, " tissues.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " inputted tissue(s) have/has been in DrugComb.\n",
          ifelse(is.null(nrow(new)), 0, nrow(new)),
          " inputted tissue(s) are/is not in DrugComb.")
  return(list(n = n, old = old, new = new))
}

#' Check whether input diseases are archived in DrugComb
#'
#' Function \code{CheckDisease} compares disease IDs passed to \code{test}
#' argument with diseases archived in DrugComb (extracted to
#' \code{\link{drugcomb}} R data). It will return the total number of diseases
#' have been archived in DrugComb as well as which diseases are already or not
#' yet archived in DrugComb.
#'
#' @param disease_df A data frame . It contains tow columns:
#'   \itemize{
#'     \item \strong{id}: the NCIt ID of the disease.
#'     \item \strong{name}: the name of the disease.
#'   }
#'
#' @return A list contains 3 elements:
#' \itemize{
#'   \item \strong{n}: It indicates that there are totally n diseases archived
#'   in DrugComb database.
#'   \item \strong{old}: It lists the names and IDs of diseases that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new}: It lists the names and IDs of diseases that
#'   \emph{are not yet} in DrugComb.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
CheckDisease <- function(disease_df) {
  message("Checking  diseases...")
  if (!all(c("id", "dname") %in% colnames(disease_df))){
    stop("The input data frame must contain columns: 'id' and 'dname'")
  }
  exist <- drugcomb$disease
  n <- nrow(exist)
  old <- unique(stats::na.omit(exist[match(disease_df$id, exist$id), ]))
  new <- unique(stats::na.omit(disease_df[!disease_df$id %in% exist$id, ]))
  message("DrugComb has archived ", n, " diseases.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " inputted disease(s) have/has been in DrugComb.\n",
          ifelse(is.null(nrow(new)), 0, nrow(new)),
          " inputted disease(s) are/is not in DrugComb.")
  return(list(n = n, old = old, new = new))
}

#' Check whether input drugs are archived in DrugComb
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
#'   \item \strong{n}: It indicates that there are totally n drugs archived in
#'   DrugComb database.
#'   \item \strong{old}: It lists the CIDs of drugs that
#'   \emph{have been archived} in DrugComb.
#'   \item \strong{new}: It lists the CIDs of drugs that \emph{are not yet}
#'   in DrugComb.
#' }
#'
#' @author
#' Jing Tang \email{jing.tang@helsinki.fi}
#' Shuyu Zheng \email{shuyu.zheng@helsinki.fi}
#'
#' @export
CheckDrug <- function(cids) {
  message("Checking drugs...")
  exist <- drugcomb$drug
  n <- nrow(exist)
  m <- max(exist$id)
  test <- stats::na.omit(cids)
  old <- stats::na.omit(exist[match(cids, exist$cid), ])
  new <- cids[!cids %in% exist$cid]
  message("DrugComb has archived ", n, " drugs.\n",
          ifelse(is.null(nrow(old)), 0, nrow(old)),
          " inputted drug(s) have/has been in DrugComb.\n",
          ifelse(is.null(length(new)), 0, length(new)),
          " inputted drug(s) are/is not in DrugComb.")
  return(list(n = n, old = old, new = new, max_id = m))
}