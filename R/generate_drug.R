#' Annotate drug
#'
#' This function is a wrapper for multiple query functions. It will query
#' PubChem, ChEMBL, UniChem, and DTC to fetch information for drugs
#'
#' @param drug_names A vector of characters containing the names for the drugs
#'
#' @return A list with two data frames:
#'   \itemize{
#'     \item \strong{drug} It contain the basic information from PubChem,
#'       ChEMBL, and UniChem database.
#'     \item \strong{target} It contains the target information for durgs from
#'       DTC database.
#'   }
#'
#' @export
#'
#' @examples
#' AnnotateDrug(c("aspirin", ""))
AnnotateDrug <- function(drug_names){
  drug <- GetCid(drug_names, type = "name")
  # 1. PubChem
  # clinical trail phase
  tmp <- GetPubPhase(unique(na.omit(drug$cid)))
  drug <- dplyr::left_join(drug, tmp, by = "cid")
  # drug property
  tmp <- GetPubchemPro(unique(na.omit(drug$cid)))
  drug <- dplyr::left_join(drug, tmp, by = "cid")
  # 2. ChEMBL
  tmp <- GetChembl(unique(na.omit(drug$inchikey)))
  drug <- dplyr::left_join(drug, tmp, by = c("inchikey" = "input_id"))
  # 3. UniChem
  tmp <- GetIds(unique(na.omit(drug$inchikey)))
  drug <- dplyr::left_join(drug, tmp, by = "inchikey")
  ## Clean chembl_id
  drug$chembl_id <- apply(drug[, c("chembl_id.x", "chembl_id.y")], 1,
                          function(x){
                            tmp <- stats::na.omit(unique(x))
                            tmp <- tmp[which(tmp != "")]
                            paste(stats::na.omit(unique(x)), collapse = "; ")
                          }
  )
  drug$chembl_id[drug$chembl_id == ""] <- NA
  # Merge clinical phase
  drug$clinical_phase <- apply(drug[, grepl("phase", colnames(drug), fixed = TRUE)],
                               1,
                               function(x) {
                                 if (all(is.na(x))){
                                   return(NA)
                                 } else {
                                   return(max(x, na.rm = TRUE))
                                 }
                                 })

  drug <- drug[, c("input_id", "cid", "chembl_id", "uni_drugbank", "uni_kegg_c",
                   "inchikey", "smiles", "molecular_formula", "clinical_phase")]
  colnames(drug) <- c("input_name", "cid", "chembl_id", "drugbank_id", "kegg_c_id",
                      "inchikey", "smiles", "molecular_formula", "clinical_phase")
  # DTC for target
  chembl_id <- unique(na.omit(drug$chembl_id))
  if (length(chembl_id) > 0) {
    target <- GetTarget(chembl_id)
  } else {
    target <- NULL
    warning("There is no CHEMBL ID for the drugs. Can not query drug target information.")
  }
  target <- drug %>%
    dplyr::select(input_id, chembl_id) %>%
    unique() %>%
    dplyr::left_join(target, by = "chembl_id")
  return(list(drug_annotation = drug, target = target))
}