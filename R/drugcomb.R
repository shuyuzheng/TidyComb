################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng<shuyu.zheng@helsinki.fi>, November 2020
################################################################################

#' Identifiers in DrugComb
#'
#' A list containing the identifiers of drugs, diseases, tissues, cell_lines,
#' studies archived in DrugComb.
#'
#' This is a list of identifiers of drugs, diseases, tissues, cell_lines,
#' studies archived in DrugComb. Used to check the records currently exist in
#' DrugComb and avoid generate duplicated entries.
#'
#' The data is request from DrugComb database, using MySQL commands:
#'
#' \code{SELECT cellosaurus_accession,id FROM cell_line
#' SELECT id,cid FROM drug
#' SELECT * FROM disease
#' SELECT * FROM tissue
#' SELECT id FROM study;}
#'
#'
#' @format A list with 5 data frames:
#' \describe{
#'   \item{cell_line}{Contains 2 variables: cellosaurus_accession, id.}
#'   \item{drug}{Contains 2 variables: id, cid.}
#'   \item{disease}{Contains 2 variables: id, name.}
#'   \item{tissue}{Contains 2 variables: name, id.}
#'   \item{study}{Contains 1 variable: id.}
#' }
#'
#' @source url{https://drugcomb.fimm.fi/}
"drugcomb"