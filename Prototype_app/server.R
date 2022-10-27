source("global.R")

# Start up celery
reticulate::use_virtualenv("/venv")
clry <- reticulate::import('celery')
celery_app <- clry$Celery('app', broker=redis_url, backend=redis_url)

server <- function(input, output, session) {
  #set to 0 for running local, 1 for running in docker
  # Sys.setenv("DEMO_VERSION" = "0")
  Sys.setenv("DEMO_VERSION" = "1")
  TheTable <- reactiveValues(ID = NULL, Unfiltered = NULL, Filtered = NULL, Job = NULL)
  
  # Reading in user uploaded file
  readUserData <- reactive({
    file <- input$inFile
    ext <- tools::file_ext(file$datapath)
    if(ext == "csv"){
      userData <- read.csv(file$datapath)
    } else {
      userData <- readRDS(file$datapath)
    }
    userData
  })
  
  observeEvent(input$pcsfSubmit, {
      isolate(userData <- readUserData())
      isolate(Protein <- input$Protein)
      isolate(logFC <- input$logFC)
      isolate(P_val <- input$P_val)
      isolate(P_val_cut <- input$P_val_cut)
      isolate(species <- input$species)
      
      #converting given column names to known ones
      colnames(userData)[which(names(userData) == Protein)] <- "Protein_Identifier"
      colnames(userData)[which(names(userData) == logFC)] <- "Fold_change"
      colnames(userData)[which(names(userData) == P_val)] <- "P_value"
      
      #filtering by provided p-value
      filtUserData <- filter(userData, P_value <= input$P_val_cut)
      
      #Adding absolute value of log fold change column
      filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
      filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
      TheTable$ID <- put_data(con, filtUserData)
      print(celery_app)
      print(TheTable)
      print(TheTable$ID)
      TheTable$Job <- celery_app$send_task("filterFun",
                       kwargs = list(
                         username = Sys.getenv("SHINYPROXY_USERNAME"),
                         id = TheTable$ID, 
                         species = species
                       ))
  })

  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
      print(TheTable$Job)
      query <- parseQueryString(session$clientData$url_search)
      message(paste0("QUERY:", query))
      # if (!is.null(query))
      # subnet <- get_data(con, query)
      # visIgraph(subnet) %>% visNodes(value = 45, size = 15, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
      #   visInteraction(navigationButtons = T) %>% visLegend(zoom = F, ncol = 2, addNodes = list(
      #     list(label = "Terminal", shape = "dot"),
      #     list(label = "Steiner", shape = "triangle"))) %>%
      #   visEvents(select = "function(nodes) {
      #           Shiny.onInputChange('pcsf_node', nodes.nodes);
      #           ;}") %>%
      #   visExport(type = "png", name = "exported-network", float = "right",
      #             label = "Export PNG", background = "white", style= "")
    }})
  observeEvent(input$checkStatus, {
    # TheTable$Job <- sendData()
    if (!is.null(TheTable$Job)) {
      sendModalAlert(TheTable$Job$info)
    } else {
      sendModalAlert("No jobs are currently running.")
    }
  })
  
}