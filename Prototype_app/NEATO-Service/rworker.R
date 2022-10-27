#!/usr/bin/env Rscript  
library(rworker)
library(reticulate)

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
filterFun <- function(username,
                      id,
                      species) {
  
  # Load libraries
  library(mapDataAccess)
  library(dplyr)
  library(tidyr)
  library(mapDataAccess)
  library(visNetwork)
  library(igraph)
  library(PCSF)
  source("enrichment.R")
  source("String.R")
  
  # Set status message
  task_progress("Pulling data")
  
  # Connect to minio
  miniocon = map_data_connection("minio_config.yml")
  
  # Add shiny proxy username variable to global environment
  Sys.setenv("SHINYPROXY_USERNAME" = username)
  Sys.setenv("DEMO_VERSION" = "1")
  
  # Pull data
  filtUserData <- get_data(miniocon, id)
  message(id)
  
  # # Get tags
  # tags <- get_tags(miniocon, id)
  
  #Set staus message
  task_progress("Making Mapping Databasbe")
  mapping_db <- make_mapping_db(species)
  
  # Set status message
  task_progress("Mapping Protein Names to String IDs")
  #adding a STRING_id column (to find protein interactions) and a genename column (to find enrichment)
  filtUserData <- map_to_stringids(filtUserData, mapping_db)
  task_progress("Mapping String IDs to readable gene names")
  filtUserData$GeneName <- convert_to_genename(filtUserData$STRING_id, mapping_db)
  
  filt_id <- put_data(miniocon, filtUserData)
  
  # # Set tags
  # set_tags(miniocon, id, list("data" = tags$Dataset))
  
  # Return status
  task_progress(paste0("Load filtered data with http://localhost:8080/?data=", filt_id))
  Sys.sleep(60)
}

# Register the task with redis
consumer$task(filterFun, name = "filterFun")

# Set consumer endpoint
consumer$consume()


