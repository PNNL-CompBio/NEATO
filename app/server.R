source("global.R")

# Start up celery
reticulate::use_virtualenv("/venv")
clry <- reticulate::import('celery')
celery_app <- clry$Celery('app', broker=redis_url, backend=redis_url)

server <- function(input, output, session) {
  # Sys.setenv("DEMO_VERSION" = "0")
  # Sys.setenv("DEMO_VERSION" = "1")
  # Sys.setenv("SHINYPROXY_USERNAME" = "loganlewis51@pnnl.gov")
  TheTable <- reactiveValues(ID = NULL, Unfiltered = NULL, Filtered = NULL, Job = NULL)

  #Help Buttons
  observeEvent(input$fileHelp, {
    showModal(modalDialog(
      title = "Upload your data",
      "Take pmart midpoint file with proteomics data, or a .csv or .tsv with protein identifiers, log-fold changes, and a p-values columns.",
      easyClose = T,
      footer = NULL
    ))
  })

  observeEvent(input$testDataHelp, {
    showModal(modalDialog(
      title = "Upload your data",
      "Loads in a test dataset containing Yeast proteins. Learn more about the data in the \"About\" page",
      easyClose = T,
      footer = NULL
    ))
  })

  # Reading in user uploaded file
  readUserData <- reactive({

    validate(
      need(!is.null(input$inFile) | input$testDataButton, "No Dataset to process")
    )

    if(!is.null(input$inFile) | input$testDataButton){
    if(input$testDataButton){
      userData <- read.csv("test_data/test_yeast.csv")
    }
    else {
      file <- input$inFile
      ext <- tools::file_ext(file$datapath)
    if(ext == "csv"){
      userData <- read.csv(file$datapath)
    } else {
      userData <- readRDS(file$datapath)
    }
    }
    userData
    }
  })

  # Displays datatable of uploaded file
  toListen <- reactive({
    list(input$testDataButton,input$inFile)
  })

  observeEvent(toListen(), {
    if(!is.null(input$inFile) | input$testDataButton){
    updateSelectInput(
      session,
      "Protein",
      choices=names(readUserData()))
  }}, ignoreInit = T)

  observeEvent(toListen(), {
    if(!is.null(input$inFile) | input$testDataButton){
    updateSelectInput(
      session,
      "logFC",
      choices=names(readUserData()))
  }}, ignoreInit = T)

  observeEvent(toListen(), {
    if(!is.null(input$inFile) | input$testDataButton){
    updateSelectInput(
      session,
      "P_val",
      choices=names(readUserData()))
  }}, ignoreInit = T)

  observeEvent(toListen(), ignoreInit = T, {
    output$previewTable <- renderDataTable({
      DT::datatable(readUserData(), options = list(scrollX = TRUE))
    }, options = list(pageLength = 10))})

  # observeEvent(input$sendService) {
  #   isolate(algorithm <- input$algorithm)
  #   output$algoTesting <- renderText(paste0("The algorithm is: ", algorithm))
  # }

  observeEvent(input$sendService, {
      output$Loading <- renderText("Submitting Data...")
    }, ignoreInit = T
  )

  # reactive({
  #   if(input$sendService > 0) {
  #   output$Loading <- renderText("Submitting Data...")
  #   }
  # })

  observeEvent(input$sendService, {
      Sys.sleep(6)
      isolate(userData <- readUserData())
      isolate(Protein <- input$Protein)
      isolate(logFC <- input$logFC)
      isolate(P_val <- input$P_val)
      isolate(P_val_cut <- input$P_val_cut)
      isolate(species <- input$species)
      isolate(includedProts <- input$includedProts)
      isolate(scoreThresh <- input$scoreThresh)
      isolate(algorithm <- input$algorithm)

      validate(
        need(length(unique(c(Protein, logFC, P_val))) == 3, "Selected column names must be unique.")
      )

      validate(
        need(!("" %in% c(Protein, logFC, P_val, P_val_cut, species, includedProts, scoreThresh, algorithm)), "A parameter was left blank! Please check your parameters.")
      )

      #converting given column names to known ones
      colnames(userData)[which(names(userData) == Protein)] <- "Protein_Identifier"
      colnames(userData)[which(names(userData) == logFC)] <- "Fold_change"
      colnames(userData)[which(names(userData) == P_val)] <- "P_value"

      # output$algoTesting <- renderText(paste0("The algorithm is: ", algorithm))
      #Submiting task
      TheTable$ID <- put_data(con, userData)
      if(algorithm == "explicitGraph") {
      TheTable$Job <- celery_app$send_task("explicitGraph",
                       kwargs = list(
                         username = Sys.getenv("SHINYPROXY_USERNAME"),
                         id = TheTable$ID,
                         species = species,
                         includedProts = includedProts,
                         scoreThresh = scoreThresh
                       ))
      } else if(algorithm == "inferredGraph") {
        TheTable$Job <- celery_app$send_task("inferredGraph",
                                             kwargs = list(
                                               username = Sys.getenv("SHINYPROXY_USERNAME"),
                                               id = TheTable$ID,
                                               species = species,
                                               includedProts = includedProts,
                                               scoreThresh = scoreThresh
                                             ))
      }
      output$Submitted <- renderText("Data submitted! Click \"Check Status\" to see the status of the job.")
  })
  enrich_db <- reactive({
    if(input$pcsfSubmit >= 1){
      query <- parseQueryString(session$clientData$url_search)
      data <- get_data(con, query)
      species <- data[[3]]
      make_enrich_db(species)
    }
  })

  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
      isolate(algorithm <- input$algorithm)
      query <- parseQueryString(session$clientData$url_search)
      print(query)
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
