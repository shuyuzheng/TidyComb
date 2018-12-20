
#' Generate drug table and drug id index table
#'
#' @param cids A vector of characters contains interested drugs
#'
#' @return A list contains 2 data frame: one contains new drug's information,
#'  another one is drug id cid index table.
#' @export
#'
#' @importFrom magrittr %>%
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

  drugbank <- utils::read.csv(system.file("extdata", "drugbank.csv",
                                   package = "TidyComb"),
                              stringsAsFactors = FALSE) %>%
     dplyr::select(db_drugbank = "DrugBank.ID",
                   db_kegg_c = "KEGG.Compound.ID",
                   db_kegg_d = "KEGG.Drug.ID",
                   cid = "PubChem.Compound.ID")

  # 2. Get infomation from databases

  drug <- GetPubNames(check$new) %>%
     dplyr::left_join(GetPubchemPro(check$new), by = "cid") %>%
     dplyr::left_join(GetPubPhase(check$new), by = "cid")

  drug <- drug %>%
     dplyr::left_join(GetChembl(drug$inchikey), by = "inchikey") %>%
     dplyr::mutate(clinical_phase = max(phase, chembl_phase)) %>%
     dplyr::left_join(GetIds(drug$inchikey), by = "inchikey") %>%
     dplyr::left_join(drugbank, by = "cid")

  # 3. clean IDs
  message("Combining drug information...")
  # KEGG Coumpound ID
  drug$db_kegg_c <- as.character(drug$db_kegg_c)
  drug$uni_kegg_c[is.na(drug$uni_kegg_c)] <- ""
  drug$db_kegg_c[is.na(drug$db_kegg_c)] <- ""
  drug$db_drugbank[is.na(drug$db_drugbank)] <- ""
  drug$uni_drugbank[is.na(drug$uni_drugbank)] <- ""

  n <- nrow(drug)
  for (i in 1:n) {
    if (drug$db_kegg_c[i] != drug$uni_kegg_c[i]) {
      ids <- c(drug$db_kegg_c[i], drug$uni_kegg_c[i])
      drug$kegg_c[i] <- paste(ids[ids != ""], collapse = "; ")
    } else {
      drug$kegg_c[i] <- drug$db_kegg_c[i]
    }
    if (drug$uni_drugbank[i] != drug$db_drugbank[i]) {
      ids <- c(drug$uni_drugbank[i], drug$db_drugbank[i])
      drug$drugbank[i] <- paste(ids[ids != ""], collapse = "; ")
    } else {
      drug$drugbank[i] <- drug$db_drugbank[i]
    }
  }

  drug$drugbank[drug$drugbank == ""] <- NA
  drug$kegg_c[drug$kegg_c == ""] <- NA

  message("Generating drug ID...")
  drug$id <- seq(check$n + 1, length.out = n)
  drug <- dplyr::select(drug, name, id, inchikey, smiles, cid, chembl_id,
                 drugbank_id = "drugbank", kegg_c, kegg_d = "db_kegg_d",
                 molecular_formula, clinical_phase)

  # 7. Generate drug ID index table.
  message("Generating 'drug_id' index table...")
  drug.id <- dplyr::select(drug, id, cid) %>%
     rbind.data.frame(check$old)

  return(list(drug = drug, drug_id = drug.id))
 }
