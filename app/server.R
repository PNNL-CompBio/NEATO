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
  
  # Displays datatable of uploaded file
  observeEvent(input$inFile, {
    output$previewTable <- renderDataTable({
      userData <- readUserData()
      DT::datatable(userData, options = list(scrollX = TRUE))
    }, options = list(pageLength = 10))})
  
  observeEvent(input$sendService, {
      isolate(userData <- readUserData())
      isolate(Protein <- input$Protein)
      isolate(logFC <- input$logFC)
      isolate(P_val <- input$P_val)
      isolate(P_val_cut <- input$P_val_cut)
      isolate(species <- input$species)
      isolate(includedProts <- input$includedProts)
      isolate(scoreThresh <- input$scoreThresh)
      isolate(algorithm <- input$algorithm)
      
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
      TheTable$Job <- celery_app$send_task(algorithm,
                       kwargs = list(
                         username = Sys.getenv("SHINYPROXY_USERNAME"),
                         id = TheTable$ID,
                         species = species,
                         includedProts = includedProts,
                         scoreThresh = scoreThresh
                       ))
  })
  enrich_db <- reactive({
    if(input$pcsfSubmit >= 1){
      isolate(species <- input$species)
      make_enrich_db(species)
    }
  })

  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
      isolate(algorithm <- input$algorithm)
      query <- parseQueryString(session$clientData$url_search)
      data <- get_data(con, query)
      subnet <- data[[2]]
      if(algorithm == "inferredGraph") {
      visIgraph(subnet) %>% visNodes(value = 45, size = 15, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
        visInteraction(navigationButtons = T) %>% visLegend(zoom = F, ncol = 2, addNodes = list(
          list(label = "Terminal", shape = "dot"),
          list(label = "Steiner", shape = "triangle"))) %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('pcsf_node', nodes.nodes);
                ;}") %>%
        visExport(type = "png", name = "exported-network", float = "right",
                  label = "Export PNG", background = "white", style= "")
      } else {
        visIgraph(subnet) %>% visNodes(value = 45, size = 15, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
          visInteraction(navigationButtons = T) %>% visLegend(zoom = F, ncol = 2) %>%
          visEvents(select = "function(nodes) {
                Shiny.onInputChange('pcsf_node', nodes.nodes);
                ;}") %>%
          visExport(type = "png", name = "exported-network", float = "right",
                    label = "Export PNG", background = "white", style= "")
      }
    }})
  observeEvent(input$checkStatus, {
    # TheTable$Job <- sendData()
    if (!is.null(TheTable$Job)) {
      sendModalAlert(TheTable$Job$info)
    } else {
      sendModalAlert("No jobs are currently running.")
    }
  })
  
  
  observeEvent(input$pcsf_node,{
    #Selecing a cluster based on click
    query <- parseQueryString(session$clientData$url_search)
    data <- get_data(con, query)
    subnet <- data[[2]]
    pcsfSelGroup <- V(subnet)$group[V(subnet)$name == input$pcsf_node]
    selNodeList <- V(subnet)$name[V(subnet)$group == pcsfSelGroup]
    visNetworkProxy("pcsf") %>% visSelectNodes(selNodeList)
  })
  
  enrich_table_pcsf <- reactive({
    isolate(enrich_db <- enrich_db())
    query <- parseQueryString(session$clientData$url_search)
    data <- get_data(con, query)
    subnet <- data[[2]]
    filtUserData <- data[[1]]
    
    #making selection based on button click
    if(input$bgChoicePCSF == "full"){
      background <- unique(c(filtUserData$GeneName, V(subnet)$name))
    }
    else {
      background <- V(subnet)$name
    }
    #filters enrichment table based on selected cluster
    if (!is.null(input$pcsf_node)) {
      pcsfSelGroup <- V(subnet)$group[V(subnet)$name == input$pcsf_node]
      clust_enrich = find_enrich(background, V(subnet)$name[which(V(subnet)$group == pcsfSelGroup)], enrich_db)
    }
    else {
      clust_enrich = find_enrich(background, V(subnet)$name, enrich_db)
    }
    clust_enrich
  })
  
  #edits table for display clarity
  edited_table_pcsf <- reactive({
    if(input$pcsfSubmit >= 1) {
      enrich_table <- enrich_table_pcsf()
      edited_table <- edit_e_table(enrich_table, input$pvalFilterPCSF, input$BHpvalFilterPCSF)
      edited_table
    }
  })
  
  #displays table
  output$pcsf_enrich <- renderDataTable({
    edited_table_pcsf <- edited_table_pcsf()
  }, options = list(pageLength = 10), rownames = F)
  
  observeEvent(input$pcsfSubmit, {
    #creating enrichment table download button
    output$downPCSFButtonFull <- renderUI({downloadButton("downloadPCSFFull", "Download Full Table")})
    output$downPCSFButtonCurrent <- renderUI({downloadButton("downloadPCSFCurrent", "Download Current Table")})
  })
  
  #operates full or current table Download buttons
  output$downloadPCSFFull <- downloadHandler(
    filename = "PCSFEnrichmentTableFull.csv",
    content = function(file) {
      write.csv(enrich_table_pcsf(), file, row.names = F)
    }
  )
  output$downloadPCSFCurrent <- downloadHandler(
    filename = "PCSFEnrichmentTableSection.csv",
    content = function(file) {
      write.csv(edited_table_pcsf(), file, row.names = F)
    }
  )
  
}