# TidyComb
# Function for generate "drug" table which will be uploaded to DrugComb
# Copyright: Shuyu Zheng

#' Generate drug table and drug id index table
#'
#' @param cids A vector of characters contains interested drugs
#'
#' @return A list contains 2 data frame: one contains new drug's information,
#'  another one is drug id cid index table.
#' @export
#'
#' @examples
#' input <- utils::read.csv(system.file("extdata",
#'                               "template.csv",
#'                                package = "TidyComb"),
#'                   stringsAsFactors = FALSE)
#' cids <- na.omit(unique(c(input$drug_row_cid, input$drug_col_cid)))
#' drug <- GenerateDrug(cids)
GenerateDrug <- function(cids) {
  # 1. Check whether drug have been archived in DrugComb database.
  check <- CheckDrug(cids = cids)

  # 2. Get infomation from databases
  # PubChem
  drug <- dplyr::left_join(GetPubNames(check$new),
                           GetPubchemPro(check$new),
                           by = "cid")
  drug <- dplyr::left_join(drug, GetPubPhase(check$new), by = "cid")

  # UniChem
  drug <- dplyr::left_join(drug, GetIds(stats::na.omit(drug$inchikey)),
                           by = "inchikey")

  # ChEMBL
  chembl_id <- stats::na.omit(drug$chembl_id)
  if (length(chembl_id) != 0) {
  drug <- dplyr::left_join(drug, GetChemblPhase(chembl_id),
                           by = "chembl_id")
  drug$chembl_phase[is.na(drug$chembl_phase)] <- 0
  drug$clinical_phase <- apply(drug[, grepl(".*phase$", colnames(drug))],
                               1, max)
  } else {
    drug$clinical_phase <- drug$phase
  }

  # DrugBank
  drugbank <- utils::read.csv(system.file("extdata", "drugbank.csv",
                                          package = "TidyComb"),
                              stringsAsFactors = FALSE)
  drugbank <- dplyr::select(drugbank, db_drugbank = "DrugBank.ID",
                            db_kegg_c = "KEGG.Compound.ID",
                            db_kegg_d = "KEGG.Drug.ID",
                            cid = "PubChem.Compound.ID")
  drug <-  dplyr::left_join(drug, drugbank, by = "cid")

  # 3. clean IDs
  message("Combining drug information...")
  # KEGG Coumpound ID
  # drug$db_kegg_c <- as.character(drug$db_kegg_c)
  # drug$uni_kegg_c[is.na(drug$uni_kegg_c)] <- ""
  # drug$db_kegg_c[is.na(drug$db_kegg_c)] <- ""
  # drug$db_drugbank[is.na(drug$db_drugbank)] <- ""
  # drug$uni_drugbank[is.na(drug$uni_drugbank)] <- ""

  drug$kegg_c <- apply(drug[, grepl(".*_kegg_c$", colnames(drug))], 1,
                       function(x){
                         paste(stats::na.omit(unique(x)), collapse = "; ")
                       }
  )


  drug$drugbank <- apply(drug[, grepl(".*_drugbank$", colnames(drug))], 1,
                         function(x){
                           paste(stats::na.omit(unique(x)), collapse = "; ")
                         }
  )

  # n <- nrow(drug)
  # for (i in 1:n) {
  #   if (drug$db_kegg_c[i] != drug$uni_kegg_c[i]) {
  #     ids <- c(drug$db_kegg_c[i], drug$uni_kegg_c[i])
  #     drug$kegg_c[i] <- paste(ids[ids != ""], collapse = "; ")
  #   } else {
  #     drug$kegg_c[i] <- drug$db_kegg_c[i]
  #   }
  #   if (drug$uni_drugbank[i] != drug$db_drugbank[i]) {
  #     ids <- c(drug$uni_drugbank[i], drug$db_drugbank[i])
  #     drug$drugbank[i] <- paste(ids[ids != ""], collapse = "; ")
  #   } else {
  #     drug$drugbank[i] <- drug$db_drugbank[i]
  #   }
  # }
  #
  drug$drugbank[drug$drugbank == ""] <- NA
  drug$kegg_c[drug$kegg_c == ""] <- NA
  drug$db_kegg_d[drug$db_kegg_d == ""] <- NA

  message("Generating drug ID...")
  drug$id <- seq(check$n + 1, length.out = nrow(drug))
  drug <- dplyr::select(drug, name, id, inchikey, smiles, cid, chembl_id,
                        drugbank_id = "drugbank", kegg_c, kegg_d = "db_kegg_d",
                        molecular_formula, clinical_phase)

  # 7. Generate drug ID index table.
  message("Generating 'drug_id' index table...")
  drug_id <- dplyr::select(drug, id, cid)
  drug_id <- rbind.data.frame(drug_id, check$old)

  return(list(drug = drug, drug_id = drug_id))
}
