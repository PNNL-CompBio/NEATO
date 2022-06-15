library(leapR)

make_enrich_db <- function(species){
  if(species == 9606){ #hard coded for now, will change when enrichment database files are more structured and final
    enrich_db <- read_gene_sets("Enrichment_db/GO_Biological_Process_2021.txt")
  } else {
    enrich_db <- read_gene_sets("Enrichment_db/Mouse_Gene_Atlas_fixed.txt")
  }
  return(enrich_db)
}



find_enrich <- function(background, targets, species){
  enrich_db <- make_enrich_db(species)
  leapR(geneset=enrich_db, enrichment_method="enrichment_in_sets",
                     background=background,
                     targets=targets)}