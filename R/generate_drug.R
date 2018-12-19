
#' Generate drug table and drug id index table
#'
#' @param cids A vector of characters contains interested drugs
#'
#' @return A list contains 2 data frame: one contains new drug's information,
#'  another one is drug id cid index table.
#' @export
#'
#' @examples
#' input <- read.csv(system.file("extdata",
#'                               "template.csv",
#'                                package = "TidyComb"),
#'                   stringsAsFactors = FALSE)
#' cids <- unique(c(input$drug_row_cid, input$drug_col_cid)) %>% na.omit()
#' drug <- GenerateDrug(cids)
GenerateDrug <- function(cids) {
  # 1. Check whether drug have been archived in DrugComb database.
  check <- CheckDrug(cids = cids)

  # 2. Get name, synonyms, and properties from PubChem
  pub.info <- GetPubNames(check$new) %>%
     dplyr::left_join(GetPubchemPro(check$new), by = "cid") %>%
     dplyr::left_join(GetPubPhase(check$new), by="cid")

  # 3. Get Chembl ID, clinical phase from ChEMBL.
  pub.chembl <- GetChembl(pub.info$inchikey, quiet = TRUE) %>%
     dplyr::right_join(pub.info, by = "inchikey") %>%
     dplyr::mutate(clinical_phase = max(phase, chembl_phase))

  # 4. Get DrugBank ID, KEGG ID from UniChem.
  pub.chembl.uni <- GetIds(pub.info$inchikey) %>%
    dplyr::right_join(pub.chembl, by = "inchikey")

  # 5. Get DrugBank ID, KEGG ID from DrugBank.
  drugbank <- read.csv(system.file("extdata",
                                   "drugbank.csv",
                                   package = "TidyComb"),
                       stringsAsFactors = FALSE)
  pub.chembl.uni.db <- drugbank %>%
     dplyr::select(db_drugbank = "DrugBank.ID",
                   db_kegg_c = "KEGG.Compound.ID",
                   db_kegg_d = "KEGG.Drug.ID",
                   cid = "PubChem.Compound.ID") %>%
     dplyr::right_join(pub.chembl.uni, by = c("cid" = "cid"))

  # 6. clean IDs
  # KEGG Coumpound ID
  kegg.compound <- dplyr::select(pub.chembl.uni.db, cid, uni_kegg_c, db_kegg_c)
  kegg.compound$db_kegg_c <- as.character(kegg.compound$db_kegg_c)
  kegg.compound$uni_kegg_c[is.na(kegg.compound$uni_kegg_c)] <- ""
  kegg.compound$db_kegg_c[is.na(kegg.compound$db_kegg_c)] <- ""

  n <- nrow(kegg.compound)
  for (i in 1:n) {
    if (kegg.compound$db_kegg_c[i] != kegg.compound$uni_kegg_c[i]) {
      ids <- c(kegg.compound$db_kegg_c[i],
               kegg.compound$uni_kegg_c[i])
      kegg.compound$kegg_c[i] <- paste(unique(ids[ids != ""]), collapse = "; ")

    } else {
      kegg.compound$kegg_c[i] <- kegg.compound$db_kegg_c[i]
    }
  }
  kegg.compound.join <- dplyr::select(kegg.compound, cid, kegg_c)

  # DrugBank IDs
  db.id <- dplyr::select(pub.chembl.uni.db, cid, uni_drugbank, db_drugbank)
  db.id$db_drugbank[is.na(db.id$db_drugbank)] <- ""
  db.id$uni_drugbank[is.na(db.id$uni_drugbank)] <- ""

  n <- nrow(db.id)
  for (i in 1:n) {
    if (db.id$uni_drugbank[i] != db.id$db_drugbank[i]) {
      ids <- c(db.id$uni_drugbank[i], db.id$db_drugbank[i])
      db.id$drugbank[i] <- paste(ids[ids != ""],
                                 collapse = "; ")
    } else {
      db.id$drugbank[i] <- db.id$db_drugbank[i]
    }
  }
  db.id.join <- dplyr::select(db.id, cid, drugbank)

  drug.combine <- pub.chembl.uni.db %>%
    dplyr::left_join(db.id.join, by = "cid") %>%
    dplyr::left_join(kegg.compound.join, by = "cid")

  drug.combine$id <- seq(check$n + 1, length.out = nrow(drug.combine))
  drug <- dplyr::select(drug.combine, name, id,
                  inchikey, smiles, cid, chembl_id,
                 drugbank_id = "drugbank", kegg_c, kegg_d = "db_kegg_d",
                 molecular_formula, clinical_phase)

  # 7. Generate drug ID index table.
  drug.id <- dplyr::select(drug, id, cid) %>%
     rbind.data.frame(check$old)

  return(list(drug = drug, drug.id = drug.id))
 }
