library(shiny)
library(shinyWidgets)
# library(shinyjs)
library(shinyBS)
library(dplyr)
library(tidyr)
# library(RCurl)
# library(STRINGdb)
library(mapDataAccess)
library(visNetwork)
library(igraph)
library(shinycssloaders)
library(shinydashboard)
library(PCSF)
library(DT)
# library(pmartRdata)
source("enrichment.R")
source("String.R")


# Register url, this is running in another docker container alongside this one
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

sendModalAlert <- function(message = "") {
  showModal(modalDialog(
    HTML(paste0('<span style="font-size: 22px;">', message, '</span>')),
    title = "", size = "s", easyClose = TRUE
  ))
}
con <- map_data_connection("minio_config.yml")





# postFormSmart <- function(uri, ..., .params = list(), .opts = curlOptions(url = uri),
#                           curl = getCurlHandle(), style = 'HTTPPOST',
#                           .encoding = integer(), binary = NA, .checkParams = TRUE,
#                           .contentEncodeFun = curlEscape){
#
#   res = postForm(uri, ..., .params = .params, .opts = .opts,
#                  curl = curl, style = style,
#                  .encoding = .encoding, binary = binary, .checkParams = .checkParams,
#                  .contentEncodeFun = .contentEncodeFun)
#
#
#   suppressWarnings( if(grepl("The document has moved", res)){
#
#     begin <- regexpr("href",res)+6
#     mys2=substr(res, begin, 10000000)
#     end <- regexpr('"',mys2)-1
#     uriNew = substr(mys2, 1, end)
#
#     res=postForm(uriNew, ..., .params = .params, .opts = .opts,
#                  curl = curl, style = style,
#                  .encoding = .encoding, binary = binary, .checkParams = .checkParams,
#                  .contentEncodeFun = .contentEncodeFun)
#   } )
#
#   return(res)
#
# }














