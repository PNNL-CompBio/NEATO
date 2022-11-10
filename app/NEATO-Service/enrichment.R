library(leapR)

make_enrich_db <- function(species){
  if(species == 9606){ #hard coded for now, will change when enrichment database files are more structured and final
    if(Sys.getenv("DEMO_VERSION") == "0"){
      enrich_db <- read_gene_sets("Enrichment_db/GO_Biological_Process_2021.txt")
    }
    else if(Sys.getenv("DEMO_VERSION") == "1"){
      enrich_db <- read_gene_sets("/edata/GO_Biological_Process_2021.txt")
    }
  } else if(species == 10090){
    if(Sys.getenv("DEMO_VERSION") == "0"){
      enrich_db <- read_gene_sets("Enrichment_db/Mouse_Gene_Atlas_fixed.txt")
    }
    else if(Sys.getenv("DEMO_VERSION") == "1"){
      enrich_db <- read_gene_sets("/edata/Mouse_Gene_Atlas_fixed.txt")
    }
  }
  return(enrich_db)
}



find_enrich <- function(background, targets, enrich_db){
  leapR(geneset=enrich_db, enrichment_method="enrichment_in_sets",
                     background=background,
                     targets=targets)}

edit_e_table <- function(enrich_table, pvalFilt, BHpvalFilt){
  enrich_table <- enrich_table[c("ingroup_n", "pvalue", "BH_pvalue")]
  enrich_table <- filter(enrich_table, ingroup_n > 0) %>% filter(pvalue <= pvalFilt) %>% filter(BH_pvalue <= BHpvalFilt)
  enrich_table <- data.frame(`Enrichment Terms` = row.names(enrich_table), enrich_table)
  colnames(enrich_table) <- c("Enrichment Terms", "Number in Ingorup", "P-value", "BH Adjusted Pvalue")
  rownames(enrich_table) <- NULL
  enrich_table <- enrich_table[order(enrich_table$`P-value`), ]
  return(enrich_table)
}