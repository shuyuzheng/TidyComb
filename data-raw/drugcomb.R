cell_line <- read.csv("data-raw/cell_line.csv",
                      stringsAsFactors = FALSE,
                      header = FALSE,
                      col.names = c("cellosaurus_accession", "id"))
drug <- read.csv("data-raw/drug.csv",
                 stringsAsFactors = FALSE,
                 header = FALSE,
                 col.names = c("id", "cid"))
disease <- read.csv("data-raw/disease.csv",
                    stringsAsFactors = FALSE,
                    header = FALSE,
                    col.names = c("name", "id"))
tissue <- read.csv("data-raw/tissue.csv",
                   stringsAsFactors = FALSE,
                   header = FALSE,
                   col.names = c("id", "name"))
study <- read.csv("data-raw/study.csv",
                  stringsAsFactors = FALSE,
                  header = FALSE,
                  col.names = "id")
drugcomb <- list(cell_line = cell_line,
                 drug = drug,
                 disease = disease,
                 tissue = tissue,
                 study = study)
usethis::use_data(drugcomb, overwrite = TRUE)
