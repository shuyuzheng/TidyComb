################################################################################
# Copyright Shuyu Zheng and Jing Tang - All Rights Reserved
# Unauthorized copying of this file, via any medium is strictly prohibited
# Proprietary and confidential
# Written by Shuyu Zheng <shuyu.zheng@helsinki.fi>, November 2020
################################################################################

#' TidyComb: Preparing drug combination data for DrugComb
#'
#' This package wrapped the functions that could be used to tidy drug
#' combination data from various sources into the format that is required by
#' uploading to DrugComb database.
#'
#' The TidyComb package provides three categories of important
#' functions :
#' \enumerate{
#'   \item Retrieving data from external databases, like PubChem, Chembl, ...;
#'   \item Estimating cell line's sensitivity to drugs. Methods includs fitting
#'   and visaulizing dose-response curve, ;
#'   \item Calculating drug combination effect surface.
#' }
#'
#' @section Retrieving data from databases:
#' The following functions are provided to
#'
#' @docType package
#' @name TidyComb
NULL