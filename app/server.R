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
  

  #reading in the appropriate database files based on species
  mapping_db <- reactive({
    isolate(species <- input$species)
    make_mapping_db(species)
  })
  
  inters_db <- reactive({
    isolate(species <- input$species)
    make_inters_db(species)
  })
  
  enrich_db <- reactive({
    isolate(species <- input$species)
    make_enrich_db(species)
  })
  
 
  
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

  #filtering user data
  filtUserData <- reactive({
    isolate(userData <- readUserData())
    isolate(Protein <- input$Protein)
    isolate(logFC <- input$logFC)
    isolate(P_val <- input$P_val)
    isolate(P_val_cut <- input$P_val_cut)
    isolate(mapping_db <- mapping_db())

    #converting given column names to known ones
    colnames(userData)[which(names(userData) == Protein)] <- "Protein_Identifier"
    colnames(userData)[which(names(userData) == logFC)] <- "Fold_change"
    colnames(userData)[which(names(userData) == P_val)] <- "P_value"
    
    #filtering by provided p-value
    filtUserData <- filter(userData, P_value <= input$P_val_cut)

    #Adding absolute value of log fold change column
    filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
    filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
    
    #adding a STRING_id column (to find protein interactions) and a genename column (to find enrichment)
    filtUserData <- map_to_stringids(filtUserData, mapping_db)
    filtUserData$GeneName <- convert_to_genename(filtUserData$STRING_id, mapping_db)

    filtUserData
  })

  # Processing data, returns data for nodes and edges
  graph <- reactive({
    if (input$submit >= 1) {
      #isolate makes it so element isn't reloaded when values change
      isolate(filtUserData <- filtUserData())
      isolate(includedProts <- as.numeric(input$includedProts))
      isolate(scoreThresh <- as.numeric(input$scoreThresh))
      isolate(inters_db <- inters_db())
      isolate(mapping_db <- mapping_db())

      #Grabbing columns for nodes
      hits <- filtUserData[1:includedProts, "STRING_id"]

      print("getting interactions")
      edges <- get_interactions(hits, scoreThresh, inters_db)
      
      #Currently proteins with no interations are not included in the final graph
      #here I could manually add them back by checking the uploaded proteins vs the ones in edges now?
      
      #converting STRING_ids to genenames
      edges$names1 <- convert_to_genename(edges$protein1, mapping_db)
      edges$names2 <- convert_to_genename(edges$protein2, mapping_db)
      
      #label columns appropriately
      colnames(edges) <- c("protein1", "protein2", "width", "from", "to")
      edges <- edges[,c("from", "to", "width", "protein1", "protein2")]
      edges$width <- as.numeric(edges$width)/100
      
      #convert to igraph object
      graph <- graph_from_data_frame(edges, directed = FALSE)
      
      #cluster with edge betweenness algorithm
      cluster <- cluster_edge_betweenness(graph)
      V(graph)$group <- cluster$membership
      return(graph)
    }
  })

#STRING NETWORK

  # Final processing and creating visNetwork on button press
  output$network <- renderVisNetwork({
    if(input$submit >= 1) {
      graph <- graph()
      visIgraph(graph) %>%
        visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
        visInteraction(navigationButtons = T) %>% visIgraphLayout() %>% visLegend(zoom = F, ncol = 2) %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('sel_node', nodes.nodes);
                ;}") %>%
        visExport(type = "png", name = "exported-network", float = "right",
                  label = "Export PNG", background = "white", style= "")
    }
  })

  # Selects cluster of clicked node
  observeEvent(input$sel_node,{
    #Selecing a cluster based on click
    graph <- graph()
    stringSelGroup <- V(graph)$group[V(graph)$name == input$sel_node]
    selNodeList <- V(graph)$name[V(graph)$group == stringSelGroup]
    visNetworkProxy("network") %>% visSelectNodes(selNodeList)
  })

#ENRICHMENT TABLE

  # Creates enrichment table
  enrich_table <- reactive({
    isolate(filtUserData <- filtUserData())
    isolate(enrich_db <- enrich_db())
    graph <- graph()
    #selects backgroud for enrichment based on user selection
    if(input$bgChoice == "full"){
      background <- filtUserData$GeneName
    }
    else {
      background <- V(graph)$name
    }
    
    #filters enrichment table based on selected cluster
    if (!is.null(input$sel_node)) {
      stringSelGroup <- V(graph)$group[V(graph)$name == input$sel_node]
      clust_enrich = find_enrich(background, V(graph)$name[which(V(graph)$group == stringSelGroup)], enrich_db)
    }
    else {
      clust_enrich = find_enrich(background, V(graph)$name, enrich_db)
  }
    clust_enrich
  })
  
  #edits table for display clarity
  edited_table <- reactive({
    if(input$submit >= 1) {
      enrich_table <- enrich_table()
      edited_table <- edit_e_table(enrich_table, input$pvalFilter, input$BHpvalFilter)
      edited_table
    }
  })
  #displays enrichment table
  output$enrich <- renderDataTable({
    edited_table <- edited_table()
  }, options = list(pageLength = 10), rownames = F)
  
  #creates download buttons
  observeEvent(input$submit, {
    #creating enrichment table download button
    output$downButtonFull <- renderUI({downloadButton("downloadFull", "Download Full Enrichment")})
    output$downButtonCurrent <- renderUI({downloadButton("downloadCurrent", "Download Current Enrichment")})
  })
  
  #downloads edited or full table based on button click
  output$downloadFull <- downloadHandler(
    filename = "EnrichmentTableFull.csv",
    content = function(file) {
      write.csv(enrich_table(), file, row.names = F)
    }
  )
  output$downloadCurrent <- downloadHandler(
    filename = "EnrichmentTableSection.csv",
    content = function(file) {
      write.csv(edited_table(), file, row.names = F)
    }
  )


#PCSF NETWORK
  subnet <- reactive({
    #checks for button click
    if (input$pcsfSubmit >= 1) {

    isolate(filtUserData <- filtUserData())
    isolate(includedProts <- input$includedProts)
    isolate(mapping_db <- mapping_db())
    isolate(inters_db <- inters_db())

    #Grabbing top log-fold change proteins
    terminals <- filtUserData$Fold_change_abs[1:includedProts]
    hits <- filtUserData[1:includedProts, ]
    names(terminals) = as.character(hits$STRING_id)

    #making network and enrichment table
    ppi <- prepare_interactome(inters_db)
    print(terminals)
    subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)

    #Clustering and making names readable
    clusters <- cluster_edge_betweenness(subnet)
    V(subnet)$group <- clusters$membership
    V(subnet)$name <- convert_to_genename(V(subnet)$name, mapping_db)
    E(subnet)$value <- E(subnet)$weight
    V(subnet)$shape <- ifelse(V(subnet)$type == "Terminal", "dot", "triangle")

    subnet
    }})

  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
    subnet <- subnet()
    visIgraph(subnet) %>% visNodes(value = 45, size = 15, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
      visInteraction(navigationButtons = T) %>% visLegend(zoom = F, ncol = 2, addNodes = list(
        list(label = "Terminal", shape = "dot"),
        list(label = "Steiner", shape = "triangle"))) %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('pcsf_node', nodes.nodes);
                ;}") %>%
      visExport(type = "png", name = "exported-network", float = "right",
                                    label = "Export PNG", background = "white", style= "")
  }})

  observeEvent(input$pcsf_node,{
    #Selecing a cluster based on click
    subnet <- subnet()
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
    
  enrich_table_pcsf <- reactive({
    isolate(filtUserData <- filtUserData())
    isolate(enrich_db <- enrich_db())
    subnet <- subnet()
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


#SPRAS TESTING
  #this section is for trying to test running spras from the R shiny app
  output$sprasNet <- renderVisNetwork({
    if(input$sprasSubmit >= 1) {
      #loading in needed inputs
      isolate(proData <- readUserData())
      isolate(species <- input$species)
      isolate(inters_db <- inters_db())
      isolate(mapping_db <- mapping_db())
      isolate(includedProts <- as.numeric(input$includedProts))
      isolate(scoreThresh <- as.numeric(input$scoreThresh))
      
      #grabbing the top absolute values log-fold change proteins, grabbing needed columns
      # hits <- filtUserData[, c("STRING_id", "Fold_change_abs")]
      proData <- map_to_stringids(proData, mapping_db)
      
      #Getting interactions from the string_db that contain 'hits' proteins
      # edges <- get_interactions(proData$STRING_id, scoreThresh, inters_db)
      colnames(inters_db) <- c("protein1", "protein2", "cost")
      inters_db$cost <- as.character((1000 - inters_db$cost)/1000)
      print(head(inters_db))
      
      #changing to appropriate column names and adding (arbitrary) sources and targets columns
      # colnames(hits) <- c("NODEID", "prize")
      # hits$sources <- rep("True", times = nrow(hits))
      # hits$targets <- rep("True", times = nrow(hits))
      outFile <- proData[,c('STRING_id','Fold_change_Healthy Control_vs_Severe')]
      colnames(outFile) <- c('NODEID', 'prize')
      outFile$sources <- rep("True", times = nrow(outFile))
      outFile$targets <- rep("True", times = nrow(outFile))
      outFile$prize <- as.character(abs(outFile$prize))
      
      print(head(outFile))
      
      #saving nodes and edges to files
      if(Sys.getenv("DEMO_VERSION") == "0"){
      write_tsv(outFile, "/Users/lewi052/enrich_proj/MAP/Prototype_app/spras/input2/shinyPrizes.tsv")
      write_tsv(inters_db, "/Users/lewi052/enrich_proj/MAP/Prototype_app/spras/input2/shinyEdges.tsv", col_names = F)
      }
      else if(Sys.getenv("DEMO_VERSION") == "1"){
      write_tsv(outFile, "/app/spras/input2/shinyPrizes.tsv")
      write_tsv(inters_db, "/app/spras/input2/shinyEdges.tsv", col_names = F)
      }
      
      #running the snakemake command
      # system("cd spras")
      if(Sys.getenv("DEMO_VERSION") == "0"){
        setwd("/Users/lewi052/enrich_proj/MAP/Prototype_app/spras")
        system("snakemake --cores 1 --configfile config/config_copy.yaml")
      }
      else if(Sys.getenv("DEMO_VERSION") == "1"){
      setwd("/app/spras")
      system("export https_proxy=\"http://proxy01.pnl.gov:3128\"")
      system("export HTTPS_PROXY=\"http://proxy01.pnl.gov:3128\"")
      system("snakemake --cores 1 --configfile /app/spras/config/config_copy.yaml")
      # system("docker run reedcompbio/omics-integrator-2 OmicsIntegrator -e /app/spras/input2/shinyEdges.tsv -p /app/spras/input2/shinyEdges.tsv -o /app/spras/output2")
      }
      
      #reading in the snakemake output and converting to igraph graph object
      if(Sys.getenv("DEMO_VERSION") == "0"){
      # outData <- read.table("/Users/lewi052/enrich_proj/MAP/Prototype_app/spras/output2/data0-omicsintegrator2-params-EHHWPMD/pathway.txt")
      }
      else if(Sys.getenv("DEMO_VERSION") == "1"){
      #outData <- read.table("/app/spras/output/data0-omicsintegrator2-params-EHHWPMD/pathway.txt")
      outData <- read.table("/app/spras/output2/data0-omicsintegrator2-params-IV3IPCJ/pathway.txt")
      # outData <- read.csv("/out_data/testDataGraph/graph.tsv", sep = "\t")
      }
      graph <- graph_from_data_frame(outData)
      
      #clustering graph and changing node names to readable gene names
      clusters <- cluster_edge_betweenness(graph)
      V(graph)$group <- clusters$membership
      V(graph)$name <- convert_to_genename(V(graph)$name, mapping_db)

      #visualizing graph
      visIgraph(graph) %>% visNodes(value = 45, size = 15, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
        visInteraction(navigationButtons = T) %>% visLegend(zoom = F, ncol = 2) %>%
        visExport(type = "png", name = "exported-network", float = "right",
                  label = "Export PNG", background = "white", style= "")
    }
  })
  #stops the app when browser window is closed
  session$onSessionEnded(stopApp)
}
