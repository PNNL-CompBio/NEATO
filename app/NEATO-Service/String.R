library(readr)

make_mapping_db <- function(species) {
  if(Sys.getenv("DEMO_VERSION") == "0"){
    mapping_file <- paste0("STRING_db/", species, ".protein.aliases.v11.5.txt.gz")
  } 
  else if(Sys.getenv("DEMO_VERSION") == "1"){
    mapping_file <- paste0("/idata/", species, ".protein.aliases.v11.5.txt.gz")
  }
mapping_db <- read.csv(mapping_file, sep = "\t")
return(mapping_db)
}

make_inters_db <- function(species) {
  if(Sys.getenv("DEMO_VERSION") == "0"){
    inters_file <- paste0("STRING_db/", species, ".protein.links.v11.5.txt.gz")
  }
  else if(Sys.getenv("DEMO_VERSION") == "1"){
    inters_file <- paste0("/idata/", species, ".protein.links.v11.5.txt.gz")
  }
inters_db <- read.table(inters_file, header = T)
return(inters_db)
}

# function to map Protein Identifiers to STRING IDs

map_to_stringids <- function(filtUserData, mapping_db = mapping_db){
  mapping <- filter(mapping_db, mapping_db$alias %in% filtUserData$Protein_Identifier)
  mapping <- mapping[!duplicated(mapping$alias),]
  colnames(mapping) <- c("STRING_id","Protein_Identifier", "source")

  filtUserData <- left_join(filtUserData, mapping, by = "Protein_Identifier")

  mapped_stat <- drop_na(filtUserData)
  return(mapped_stat)
}

get_interactions <- function(hits, scoreThresh, inter_db){
  edges <- inter_db[which(inter_db$protein1 %in% hits & inter_db$protein2 %in% hits), c("protein1", "protein2", "combined_score")]
  edges <- edges[which(edges$combined_score >= scoreThresh),]
  return(edges)
}


convert_to_genename <- function(hits, map_db){
  names <- c()
  for (i in hits){
    n1 <- map_db[which(map_db$X.string_protein_id == i & map_db$source == "Ensembl_EntrezGene"), 2]
    if(length(n1) >= 1){
      names <- c(names, toupper(n1[1]))
    } else{
      names <- c(names, i)
    }
  }
  return(as.vector(names))
}

prepare_interactome <- function(inters_db = inters_db){
  interactome <- as.data.frame(inters_db[, c("protein1", "protein2", "combined_score")])
  interactome$combined_score = ((1000 - interactome$combined_score)/1000)
  interactome <- construct_interactome(interactome)
  return(interactome)
}
