GetMaxPhase <- function(cids) {
  pub.phase <- GetPubPhase(cids, quiet = FALSE)
  pub.phase$cid <- as.numeric(clinical_phase$cid)

  chemb.phase <- GetChemblPhase(cids)
  drug <- left_join(drug, clinical_phase, by = "cid")
}