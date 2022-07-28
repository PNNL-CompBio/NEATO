library(shiny)
library(shinyWidgets)
# library(shinyjs)
library(shinyBS)
library(dplyr)
library(tidyr)
# library(RCurl)
# library(STRINGdb)
library(visNetwork)
library(igraph)
library(shinycssloaders)
library(shinydashboard)
library(PCSF)
library(DT)
# library(pmartRdata)
source("enrichment.R")
source("String.R")

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














