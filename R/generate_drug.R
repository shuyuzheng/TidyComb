input <- read.csv(system.file("extdata", "template.csv", package = "TidyComb"),
stringsAsFactors = FALSE)
cids <- unique(c(input$drug_row_cid, input$drug_col_cid)) %>%
  na.omit()


GenerateDrug <- function(cids) {
 # 1. Check whether drug have been archived in DrugComb database.
 check <- CheckDrug(cids = cids)

 # 2. Get name, synonyms, and properties from PubChem
 pub.info <- GetPubNames(check$new) %>%
    dplyr::left_join(GetPubchemPro(check$new), by = "cid") %>%
    dplyr::left_join(GetPubPhase(check$new), by="cid")

chembl.info <- GetChemblPhase(pub.info$inchikey)
}
