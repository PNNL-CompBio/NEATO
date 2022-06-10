library(leapR)

enrich_db <- read_gene_sets("Enrichment_db/GO_Biological_Process_2021.txt")


find_enrich <- function(background, targets){
  leapR(geneset=enrich_db, enrichment_method="enrichment_in_sets",
                     background=background,
                     targets=targets)}