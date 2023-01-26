#!/usr/bin/env Rscript  
library(rworker)
library(reticulate)
source("String.R")

print("INFO: New Rworker thread started...")
reticulate::use_virtualenv("/venv")

if(!file.exists("redis_config.yml")){
  warning("No redis configuration found, attempting connection to default url: redis://redis1:6379")
  redis_url <- "redis://redis1:6379/0"
} else {
  redis_cfg = yaml::read_yaml("redis_config.yml")
  redis_host = redis_cfg[['host']]
  redis_url <- sprintf('redis://%s:%s/%s', 
                       redis_host, 
                       redis_cfg[['port']],
                       redis_cfg[['db']])
}

message("Setting up redis connection at:  ", redis_url)

# Instantiate Rworker object --> link between worker and task manager
consumer <- rworker(name = 'celery', workers = 1, queue = redis_url, backend = redis_url)

#' @param username The SHINYPROXY_USERNAME variable. Important for ensuring the results are returned to the user's subfolder.
#' @param id The universal unique identified (uuid) of the file.
#' @param input all input values from the NEATO shiny app. 
# Create the filter function
inferredGraph <- function(username,
                      id,
                      species,
                      includedProts,
                      scoreThresh) {
  
  # Set status message
  task_progress("Loading Packages")
  
  # Load libraries
  library(mapDataAccess)
  library(dplyr)
  library(tidyr)
  library(mapDataAccess)
  library(visNetwork)
  library(igraph)
  library(PCSF)
  source("String.R")
  
  # Set status message
  task_progress("Pulling data")
  
  # Connect to minio
  miniocon = map_data_connection("minio_config.yml")
  
  # Add shiny proxy username variable to global environment
  Sys.setenv("SHINYPROXY_USERNAME" = username)
  Sys.setenv("DEMO_VERSION" = "1")
  
  # Pull data
  userData <- get_data(miniocon, id)
  
  #filtering by provided p-value
  filtUserData <- filter(userData, P_value <= input$P_val_cut)
  
  #Adding absolute value of log fold change column
  filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
  filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
  
  # # Get tags
  # tags <- get_tags(miniocon, id)
  
  #Set staus message
  task_progress("Making Mapping Databasbe")
  mapping_db <- make_mapping_db(species)
  inters_db <- make_inters_db(species)
  
  # Set status message
  task_progress("Mapping Protein Names to String IDs")
  #adding a STRING_id column (to find protein interactions) and a genename column (to find enrichment)
  filtUserData <- map_to_stringids(filtUserData, mapping_db)
  # task_progress("Mapping String IDs to readable gene names")
  filtUserData$GeneName <- convert_to_genename(filtUserData$STRING_id, mapping_db)
  
  task_progress("Making PCSF Graph")
  terminals <- filtUserData$Fold_change_abs[1:includedProts]
  hits <- filtUserData[1:includedProts, ]
  names(terminals) = as.character(hits$STRING_id)
  
  #making PCSF network
  ppi <- prepare_interactome(inters_db)
  subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
  
  #Clustering and making names readable
  task_progress("Clustering Graph")
  clusters <- cluster_edge_betweenness(subnet)
  V(subnet)$group <- clusters$membership
  
  task_progress("Adding details to Graph")
  V(subnet)$name <- convert_to_genename(V(subnet)$name, mapping_db)
  E(subnet)$value <- E(subnet)$weight
  V(subnet)$shape <- ifelse(V(subnet)$type == "Terminal", "dot", "triangle")
  
  # # Set tags
  # set_tags(miniocon, id, list("data" = tags$Dataset))
  data_list <- list(filtUserData, subnet, species)
  
  data_id <- put_data(miniocon, data_list)
  
  # Return status
  task_progress(paste0("Load filtered data with https://map.emsl.pnnl.gov/app/neato/?data=", data_id))
  Sys.sleep(210)
}

explicitGraph <- function(username,
                          id,
                          species,
                          includedProts,
                          scoreThresh) {
  # Set status message
  task_progress("Loading Packages")
  
  library(mapDataAccess)
  library(dplyr)
  library(tidyr)
  library(mapDataAccess)
  library(visNetwork)
  library(igraph)
  library(PCSF)
  source("String.R")
  
  # Set status message
  task_progress("Pulling data")
  
  # Connect to minio
  miniocon = map_data_connection("minio_config.yml")
  
  # Add shiny proxy username variable to global environment
  Sys.setenv(SHINYPROXY_USERNAME = username)
  Sys.setenv("DEMO_VERSION" = "1")
  
  # Pull data
  userData <- get_data(miniocon, id)
  
  #filtering by provided p-value
  filtUserData <- filter(userData, P_value <= input$P_val_cut)
  
  #Adding absolute value of log fold change column
  filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
  filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
  
  # # Get tags
  # tags <- get_tags(miniocon, id)
  
  #Set status message
  task_progress("Making Mapping Databasbe")
  mapping_db <- make_mapping_db(species)
  task_progress("Making Interactions Databasbe")
  inters_db <- make_inters_db(species)
  
  # Set status message
  task_progress("Mapping Protein Names to String IDs")
  #adding a STRING_id column (to find protein interactions) and a genename column (to find enrichment)
  filtUserData <- map_to_stringids(filtUserData, mapping_db)
  filtUserData$GeneName <- convert_to_genename(filtUserData$STRING_id, mapping_db)
  
  #Grabbing columns for nodes
  hits <- filtUserData[1:includedProts, "STRING_id"]
  
  task_progress("Getting Interactions")
  edges <- get_interactions(hits, scoreThresh, inters_db)
  
  #Currently proteins with no interations are not included in the final graph
  #here I could manually add them back by checking the uploaded proteins vs the ones in edges now?
  
  #converting STRING_ids to genenames
  task_progress("Adding gene names")
  edges$names1 <- convert_to_genename(edges$protein1, mapping_db)
  edges$names2 <- convert_to_genename(edges$protein2, mapping_db)
  
  #label columns appropriately
  colnames(edges) <- c("protein1", "protein2", "width", "from", "to")
  edges <- edges[,c("from", "to", "width", "protein1", "protein2")]
  edges$width <- as.numeric(edges$width)/100
  
  #convert to igraph object
  task_progress("Creating igraph object")
  graph <- graph_from_data_frame(edges, directed = FALSE)
  
  #cluster with edge betweenness algorithm
  task_progress("Clustering graph")
  cluster <- cluster_edge_betweenness(graph)
  V(graph)$group <- cluster$membership
  
  data_list <- list(filtUserData, graph, species)
  
  data_id <- put_data(miniocon, data_list)
  testing <- get_data(miniocon, data_id)
  
  # Return status
  task_progress(paste0("Load filtered data with https://map.emsl.pnnl.gov/app/neato/?data=", data_id))
  Sys.sleep(210)
}

# Register the task with redis
consumer$task(explicitGraph, name = "explicitGraph")
consumer$task(inferredGraph, name = "inferredGraph")

# Set consumer endpoint
consumer$consume()


