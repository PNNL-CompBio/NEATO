library(readr)

mapping_db <- read.csv("/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache/9606.protein.aliases.v11.5.txt.gz", sep = "\t")

inters_db <- read.table("/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache/9606.protein.links.v11.5.txt.gz", header = T)

# function to map Protein Identifiers to STRING IDs

map_to_stringids <- function(filtUserData){
  mapping <- filter(mapping_db, mapping_db$alias %in% filtUserData$Protein_Identifier)
  mapping <- mapping[!duplicated(mapping$alias),]
  colnames(mapping) <- c("STRING_id","Protein_Identifier", "source")
  
  filtUserData <- left_join(filtUserData, mapping, by = "Protein_Identifier")
  
  mapped_stat <- drop_na(filtUserData)
  return(mapped_stat)
}

get_interactions <- function(hits, scoreThresh){
  edges <- inters_db[which(inters_db$protein1 %in% hits & inters_db$protein2 %in% hits), c("protein1", "protein2", "combined_score")]
  edges <- edges[which(edges$combined_score >= scoreThresh),]
  return(edges)
}


convert_to_hugo <- function(hits){
  names <- c()
  for (i in hits){
    n1 <- mapping_db[which(mapping_db$X.string_protein_id == i & mapping_db$source == "BioMart_HUGO"), 2]
    if(length(n1) == 1){
      names <- c(names, n1)
    } else{
      names <- c(names, i)
    }
  }
  return(as.vector(names))
}

prepare_interactome <- function(){
  interactome <- as.data.frame(inters_db[, c("protein1", "protein2", "combined_score")])
  interactome$combined_score = ((1000 - interactome$combined_score)/1000)
  interactome <- construct_interactome(interactome)
  return(interactome)
}
