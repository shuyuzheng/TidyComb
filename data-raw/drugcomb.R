cell_line <- read.csv("data-raw/cell_line.csv",
                      stringsAsFactors = FALSE)
colnames(cell_line) <- c("cellosaurus_accession", "id")
drug <- read.csv("data-raw/drug.csv",
                 stringsAsFactors = FALSE)
disease <- read.csv("data-raw/disease.csv",
                    stringsAsFactors = FALSE)
tissue <- read.csv("data-raw/tissue.csv",
                   stringsAsFactors = FALSE)
drugcomb <- list(cell_line = cell_line,
                 drug = drug,
                 disease = disease,
                 tissue = tissue)
usethis::use_data(drugcomb, overwrite = TRUE, internal = TRUE)
