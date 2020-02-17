## ----Loading_package, message=FALSE, warning=FALSE----------------------------
library("dplyr")
library("TidyComb")

## ----Input_data---------------------------------------------------------------
template <- read.csv(system.file("template.csv", package = "TidyComb"),
                     stringsAsFactors = FALSE)
head(template)

## ----Download_cellosaurus_data, message = FALSE, eval=FALSE-------------------
#  # Download "cellosaurus.html" file
#  download.file("ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml",
#                "cellosaurus.xml")

## ----Match_cell_accession, message=FALSE--------------------------------------
cells <- unique(na.omit(template$cell_line_name)) # The cell name vector must be unique and without NA
cell.match <- MatchCellAcc(cells, file = system.file("cellosaurus.xml", package = "TidyComb"))
print(cell.match)

## ----clean_cell_accession-----------------------------------------------------
cell.match.clean <- cell.match[-2, c("input_name", "cellosaurus_accession")]

## ----Generate_cell_line-------------------------------------------------------
cell <- GenerateCell(cell.match.clean$cellosaurus_accession, 
                     file = system.file("cellosaurus.xml", package = "TidyComb"))
cell$cell_line
cell_index <- cell$cell_id %>% 
  left_join(cell.match.clean, by = "cellosaurus_accession")

## -----------------------------------------------------------------------------
new_tisse <- CheckTissue(cell$cell_line$tissue)

## -----------------------------------------------------------------------------
cell$cell_line$tissue <- "liver"
new_tissue <- CheckTissue(c(cell$cell_line$tissue))

## -----------------------------------------------------------------------------
tissue_index <- new_tissue$old
tissue_index

## -----------------------------------------------------------------------------
new_disease <- CheckDisease(select(cell$cell_line, id = disease_id, dname = disease_name))
disease_index <- rbind.data.frame(new_disease$old, new_disease$new)
disease_table <- new_disease$new
disease_table

## -----------------------------------------------------------------------------
cell_line_table <- cell$cell_line %>% 
  left_join(tissue_index, by = c("tissue" = "tname")) %>% 
  rename(tissue_id = "id.y") %>% 
  select(name, synonyms, cellosaurus_accession, disease_id, id = id.x, tissue_id)
cell_line_table

## ----Get_CID------------------------------------------------------------------
drug_input <- na.omit(unique(c(template$drug_row, template$drug_col)))
drug_cid <- GetCid(drug_input, type = "name")
drug_cid

## -----------------------------------------------------------------------------
drug_cid <- drug_cid[!duplicated(drug_cid$input_id), ]
drug_cid

## -----------------------------------------------------------------------------
new_drug <- CheckDrug(drug_cid$cid)

## ----Genenrate_drug, message=FALSE--------------------------------------------
drug_name <- GetPubNames(cids = new_drug$new)
drug_info <- GetPubchemPro(cids = new_drug$new)
drug <- full_join(drug_name, drug_info)
pub_phase <- GetPubPhase(drug$cid)
drug <- left_join(drug, pub_phase, by = "cid")

## -----------------------------------------------------------------------------
chembl <- GetChembl(drug$inchikey)
drug <- left_join(drug, chembl, by = "inchikey")

## -----------------------------------------------------------------------------
drug$phase <- apply(drug[, c("phase", "chembl_phase")], 1, function(x){
  max(x, na.rm = TRUE)
  })
drug <- select(drug, -chembl_phase)

## -----------------------------------------------------------------------------
unichem <- GetIds(drug$inchikey)
drug <- left_join(drug, unichem, by = "inchikey") 
drug$chembl_id <- unlist(apply(drug[, c("chembl_id.x", "chembl_id.y")], 1,
                        function(x){
                          if (length(na.omit(x)) == 0){
                            return(NA)
                          } else{
                            return(paste(na.omit(x), collapse = "; "))
                          }
                        }))
drug <- drug[, which(!colnames(drug) %in% c("chembl_id.x", "chembl_id.y"))]

## -----------------------------------------------------------------------------
drug$id <- seq(new_drug$n + 1, length.out = nrow(drug))

## -----------------------------------------------------------------------------
drug <- select(drug, dname = "name", id, chembl_id, inchikey, smiles, cid, 
               molecular_formula, clinical_phase = "phase", #cid_m, cid_s, stitch_name, 
               drugbank_id = "uni_drugbank", kegg_id = "uni_kegg_c")
drug

## -----------------------------------------------------------------------------
drug_index <- drug_cid %>% 
  left_join(new_drug$old, by = "cid") %>% 
  left_join(select(drug, cid, id), by = "cid")
drug_index$id <- apply(drug_index[, c("id.x", "id.y")], 1, function(x){
  na.omit(x)
})
drug_index <- drug_index[, which(!colnames(drug_index) %in% c("id.x", "id.y"))]
drug_index

## ----Generate_tables, message=FALSE, warning=FALSE----------------------------
tables <- CalculateTemplate(template)
names(tables)

## -----------------------------------------------------------------------------
bad_quality <- tables$synergy$block_id[which(tables$synergy$inhibition < -200 |
                                               tables$synergy$inhibition > 200)]
bad_quality <- unique(bad_quality)
summary <- tables$summary %>% 
  mutate(quality = rep(NA, n()))
if (length(bad_quality) != 0){
  summary$quality[which(summary$block_id %in% bad_quality)] <- "bad"
}

## -----------------------------------------------------------------------------
curve <- left_join(tables$curve, drug_index, by = c("drug_row" = "input_id")) %>% 
  rename(drug_row_id = "id") %>% 
  left_join(drug_index, by = c("drug_col" = "input_id")) %>% 
  rename(drug_col_id = "id") %>% 
  select(block_id, drug_row_id, drug_col_id, b, c, d, e, model)
curve
summary <- left_join(summary, drug_index, by = c("drug_row" = "input_id")) %>% 
  rename(drug_row_id = "id") %>% 
  left_join(drug_index, by = c("drug_col" = "input_id")) %>% 
  rename(drug_col_id = "id") %>% 
  left_join(cell_index, by = c("cell_line_name" = "input_name")) %>% 
  rename(cell_line_id = "id") %>% 
  select(block_id, drug_row_id, drug_col_id, cell_line_id, conc_r_unit, 
         conc_c_unit, synergy_zip, synergy_loewe, synergy_hsa, synergy_bliss,
         ic50_row, ic50_col, ri_row, ri_col, css_row, css_col, css, S, quality)
summary

## -----------------------------------------------------------------------------
cell_line_table
tissue_index
disease_table
tables$synergy
summary
tables$surface
curve

## ----Output, eval=FALSE, message=FALSE, warning=FALSE-------------------------
#  dir.create("Upload")
#  write.csv(cell_line_table, "Upload/cell_line.csv", row.names = FALSE)
#  write.csv(tissue_index, "Upload/tissue.csv", row.names = FALSE)
#  write.csv(disease_table, "Upload/disease.csv", row.names = FALSE)
#  write.csv(drug, "Upload/drug.csv", row.names = FALSE)
#  write.csv(tables$synergy, "Upload/response.csv", row.names = FALSE)
#  write.csv(summary, "Upload/summary.csv", row.names = FALSE)
#  write.csv(tables$surface, "Upload/surface.csv", row.names = FALSE)
#  write.csv(curve, "Upload/curve.csv", row.names = FALSE)

