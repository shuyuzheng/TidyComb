GetMaxPhase <- function(cids, inchikeys) {
  pub.phase <- GetPubPhase(cids, quiet = FALSE)
  pub.phase$cid <- as.numeric(pub.phase$cid)

  chemb.phase <- GetChemblPhase(inchikeys)
  drug <- dplyr::left_join(drug, chemb.phase, by = "cid")
}