source("global.R")

server <- function(input, output, session) {
  
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
  
  filtUserData <- reactive({
    isolate(userData <- readUserData())
    isolate(Protein <- input$Protein)
    isolate(logFC <- input$logFC)
    isolate(P_val <- input$P_val)
    
    #converting given column names to known ones
    colnames(userData)[which(names(userData) == Protein)] <- "Protein_Identifier"
    colnames(userData)[which(names(userData) == logFC)] <- "Fold_change"
    colnames(userData)[which(names(userData) == P_val)] <- "P_value"
    
    filtUserData <- filter(userData, P_value <= 0.05)
    
    #Adding abolute value of log fold change column
    filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
    filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
    
    filtUserData
  })
  
  mapping_db <- reactive({
    isolate(species <- input$species)
    make_mapping_db(species)
  })
  
  inters_db <- reactive({
    isolate(species <- input$species)
    make_inters_db(species)
  })
  
  
  output$constructed <- renderText({
    if(input$construct >= 1){
    mapping_db()
    inters_db()
    "Databases have been constructed!"
  }})
  
  # Processing data, returns data for nodes and edges
  graph <- reactive({
    if (input$submit >= 1) {
      #isolate makes it so element isn't reloaded when values change
      isolate(filtUserData <- filtUserData())
      isolate(includedProts <- as.numeric(input$includedProts))
      isolate(scoreThresh <- as.numeric(input$scoreThresh))
      isolate(mapping_db <- mapping_db())
      isolate(inters_db <- inters_db())
    
      #mapping based on given identifier
      print("mapping proteins")
      mapped_stat <- map_to_stringids(filtUserData, mapping_db)

      
      # #sorting by absolute value of logFC
      # mapped_stat <- mapped_stat[order(mapped_stat$`Fold_change_abs`, decreasing = T), ]
      
      #Grabbing columns for nodes
      hits <- mapped_stat[1:includedProts, "STRING_id"]
      
      print("getting interactions")
      edges <- get_interactions(hits, scoreThresh, inters_db)
      
      edges$names1 <- convert_to_genename(edges$protein1, mapping_db)
      edges$names2 <- convert_to_genename(edges$protein2, mapping_db)
      
      colnames(edges) <- c("protein1", "protein2", "width", "from", "to")
      edges <- edges[,c("from", "to", "width", "protein1", "protein2")]
      edges$width <- as.numeric(edges$width)/100

      graph <- graph_from_data_frame(edges, directed = FALSE)
      cluster <- cluster_edge_betweenness(graph)
      V(graph)$group <- cluster$membership
      
      graph
    }
  })
  
  # Reads in User uploaded file
  observeEvent(input$inFile, {
  output$previewTable <- renderDataTable({
    userData <- readUserData()
    DT::datatable(userData, options = list(scrollX = TRUE))
    }, options = list(pageLength = 10))})

#STRING NETWORK
  
  # Final processing and creating visNetwork on button press
  output$network <- renderVisNetwork({
    if(input$submit >= 1) {
      graph <- graph()
      print("making graph")
      visIgraph(graph) %>% 
        visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
        visInteraction(navigationButtons = T) %>% visIgraphLayout(randomSeed = 123) %>% visLegend(zoom = F) %>%
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
    isolate(species <- input$species)
    graph <- graph()
    #filters enrichment table based on selected cluster
    if (!is.null(input$sel_node)) {
      stringSelGroup <- V(graph)$group[V(graph)$name == input$sel_node]
      clust_enrich = find_enrich(V(graph)$name, V(graph)$name[which(V(graph)$group == stringSelGroup)], species)
    }
    else {
      clust_enrich = find_enrich(V(graph)$name, V(graph)$name, species)
  }
    clust_enrich
  })

  #Displays enrichment table
  output$enrich <- renderDataTable({
    if(input$submit >= 1) {
    enrich_table <- enrich_table()
    # enrich_table[, c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr")]
    }
    }, options = list(pageLength = 10))

  #creates Download button after table has loaded
  observeEvent(input$submit, {
    output$downButton <- renderUI({downloadButton("downloadEnrich", "Download Table")})
  })

  #operates Download button
  output$downloadEnrich <- downloadHandler(
    filename = "EnrichmentTable.csv",
    content = function(file) {
      write.csv(enrich_table(), file, row.names = F)
    }
  )

#PCSF NETWORK
  subnet <- reactive({
    if (input$pcsfSubmit >= 1) {
    
    isolate(filtUserData <- filtUserData())
    isolate(includedProts <- input$includedProts)
    isolate(mapping_db <- mapping_db())
    isolate(inters_db <- inters_db())
    
    #Grabbing top log-fold change proteins
    terminals <- filtUserData$Fold_change_abs[1:includedProts]
    hits <- filtUserData[1:includedProts, ]
    hits <- map_to_stringids(hits, mapping_db)
    names(terminals) = as.character(hits$STRING_id)
    
    #making network and enrichment table
    ppi <- prepare_interactome(inters_db)
    length(terminals)
    subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
    
    #Clustering and making names readable
    clusters <- cluster_edge_betweenness(subnet)
    V(subnet)$group <- clusters$membership
    V(subnet)$name <- convert_to_genename(V(subnet)$name, mapping_db)
    
    subnet
    }})
  
  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
    subnet <- subnet()
    visIgraph(subnet) %>% visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
      visInteraction(navigationButtons = T) %>% visLegend(zoom = F) %>%
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
    isolate(species <- input$species)
    subnet <- subnet()
    #filters enrichment table based on selected cluster
    if (!is.null(input$pcsf_node)) {
      pcsfSelGroup <- V(subnet)$group[V(subnet)$name == input$pcsf_node]
      clust_enrich = find_enrich(V(subnet)$name, V(subnet)$name[which(V(subnet)$group == pcsfSelGroup)], species)
    }
    else {
      clust_enrich = find_enrich(V(subnet)$name, V(subnet)$name, species)
    }
    clust_enrich
  })
  
  output$pcsf_enrich <- renderDataTable({
    if(input$pcsfSubmit >= 1) {
      enrich_table <- enrich_table_pcsf()
    }
  }, options = list(pageLength = 10))
  
  observeEvent(input$pcsfSubmit, {
    #creating enrichment table download button
    output$downPCSFButton <- renderUI({downloadButton("downloadPCSF", "Download Table")})
  })
  
  #operates Download button
  output$downloadPCSF <- downloadHandler(
    filename = "PCSFEnrichmentTable.csv",
    content = function(file) {
      write.csv(enrich_table_pcsf(), file, row.names = F)
    }
  )
  #stops the app when browser window is closed
  session$onSessionEnded(stopApp)
}





