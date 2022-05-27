source("global.R")

server <- function(input, output, session) {
  
  # Creating STRINGdb database
  string_db <- reactive({
    isolate(scoreThresh <- as.numeric(input$scoreThresh))
    isolate(species <- as.numeric(input$species))
    input_directory <- paste(getwd(), "/STRINGdb_cache", sep = "")
    STRINGdb$new(version="11.5", species=species, score_threshold=scoreThresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
  })
  
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
    
    #Adding a gene symbols column to use as a label
    geneSym <- UniProt.ws::select(org.Hs.eg.db, filtUserData$Protein_Identifier, "SYMBOL","UNIPROT")
    geneSym <- geneSym[!duplicated(geneSym[,1]),]
    filtUserData <- mutate(filtUserData, geneSymbol = geneSym[,2])
    
    #Adding abolute value of log fold change column
    filtUserData <- mutate(filtUserData, Fold_change_abs = abs(filtUserData$Fold_change))
    filtUserData <- filtUserData[order(filtUserData$`Fold_change_abs`, decreasing = T), ]
    
    filtUserData
  })
  
  # Processing data, returns data for nodes and edges
  hits <- reactive({
    if (input$submit >= 1) {
      string_db <- string_db()
      
      #isolate makes it so element isn't reloaded when values change
      isolate(filtUserData <- filtUserData())
      isolate(includedProts <- as.numeric(input$includedProts))
      isolate(scoreThresh <- as.numeric(input$scoreThresh))
      isolate(species <- as.numeric(input$species))
    
      #mapping based on given identifier
      mapped_stat <- string_db$map(filtUserData, "Protein_Identifier", removeUnmappedRows = TRUE )
      
      #sorting by absolute value of logFC
      mapped_stat <- mapped_stat[order(mapped_stat$`Fold_change_abs`, decreasing = T), ]
      
      #Grabbing columns (for now grabs hard coded "Gene_Names" column, will look into best way to handle this dynamically)
      nodes <- mapped_stat[1:includedProts, c("geneSym", "STRING_id")]
      
      #(For now) Directly querying the online STIRNG database for the interactions
      identifiers = ""
      for (id in nodes$STRING_id) {
        identifiers = paste(identifiers, id, sep = "%0d")
      }
      urlStr = "https://string-db.org/api/tsv/network?"
      params <- list(identifiers = identifiers, species = species, required_score = scoreThresh)
      edges <- read_tsv(postFormSmart(urlStr, .params = params))
      
      colnames(nodes)  <- c("label", "id")
      
      # edges <- string_db$get_interactions(mapped_stat$STRING_id)
      edges <- edges[, c("stringId_A", "stringId_B", "score")]
      edges <- edges[!duplicated(edges),]
      colnames(edges) <- c("from", "to", "width")
      edges$width <- as.integer(as.numeric(edges$width)*20)

      #returns a list of the nodes and edges
      list(nodes = nodes,
           edges = edges)
    }
  })
  
  # Clusters nodes
  nodes <- reactive({
    string_db <- string_db()
    hits <- hits()
    nodes <- hits$nodes
    
    #Uses stringdb's clusters
    clusterList <- string_db$get_clusters(nodes$id)
    
    #Finds the clusters assigned to node by stringdb
    clusters <- c()
    for (x in nodes$id) {
      clusters <- c(clusters, grep(x, clusterList))
    }
    nodes$group <- clusters
    
    #Cannot have duplicated ids
    nodes <- nodes[!duplicated(nodes$id),]
    nodes <- nodes[order(nodes$group), ]
    
    #If there's only one member of the cluster, it's relabeled as "No Cluster"
    for (y in nodes$group){
      if (sum(nodes$group==y) == 1){
        nodes$group[which(nodes$group == y)] <- "No Cluster"
      }
    }
    nodes
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
      nodes <- nodes()
      hits <- hits()
      edges <- hits$edges
      
      visNetwork(nodes = nodes, edges = edges) %>% 
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
    nodes <- nodes()
    sel_node_row <- nodes[nodes$id == input$sel_node, ]
    node_cluster <- filter(nodes, group == sel_node_row$group)
    visNetworkProxy("network") %>% visSelectNodes(node_cluster$id)
  })
  
#ENRICHMENT TABLE
  
  # Creates enrichment table
  enrich_table <- reactive({
    string_db <- string_db()
    nodes <- nodes()
    #filters enrichment table based on selected cluster
    if (!is.null(input$sel_node)) {
      sel_node_row <- nodes[nodes$id == input$sel_node, ]
      nodes <- filter(nodes, group == sel_node_row$group)
    }
    enrichment <- string_db$get_enrichment(nodes$id)
    enrichment
  })
  
  #Displays enrichment table
  output$enrich <- renderDataTable({
    if(input$submit >= 1) {
    enrich_table <- enrich_table()
    enrich_table[, c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr")]
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
  res_pcsf <- reactive({
    if (input$pcsfSubmit >= 1) {
    
    isolate(filtUserData <- filtUserData())
    isolate(includedProts <- input$includedProts)
    
    #Grabbing top log-fold change proteins
    print(colnames(filtUserData))
    terminals <- filtUserData$Fold_change_abs[1:includedProts]
    names(terminals) = as.character(filtUserData$geneSymbol[1:includedProts])
    print(names(terminals))
    
    #making network and enrichment table
    ppi <- construct_interactome(STRING)
    subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
    res <- enrichment_analysis(subnet)
    }})
  
  output$pcsf <- renderVisNetwork({
    if (input$pcsfSubmit >= 1) {
    res <- res_pcsf()
    visIgraph(res$subnet) %>% visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
      visInteraction(navigationButtons = T) %>% visLegend(zoom = F) %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('pcsf_node', nodes.nodes);
                ;}") %>%
      visExport(type = "png", name = "exported-network", float = "right",
                                    label = "Export PNG", background = "white", style= "")
  }})
  
  observeEvent(input$pcsf_node,{
    #Selecing a cluster based on click
    res <- res_pcsf()
    pcsfSelGroup <<- V(res$subnet)$group[V(res$subnet)$name == input$pcsf_node]
    selNodeList <- V(res$subnet)$name[V(res$subnet)$group == pcsfSelGroup]
    visNetworkProxy("pcsf") %>% visSelectNodes(selNodeList)
  })
  
  output$pcsf_enrich <- renderDataTable({
    #displayihng enrichment table
    res <- res_pcsf()
    if (!is.null(input$pcsf_node)) {
      res$enrich <- dplyr::filter(res$enrich, res$enrich$Cluster == pcsfSelGroup)
    }
    res$enrich[, c("Cluster", "Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")]
  }, options = list(pageLength = 10))
  
  observeEvent(input$pcsfSubmit, {
    #creating enrichment table download button
    output$downPCSFButton <- renderUI({downloadButton("downloadPCSF", "Download Table")})
  })
  
  #operates Download button
  output$downloadPCSF <- downloadHandler(
    filename = "PCSFEnrichmentTable.csv",
    content = function(file) {
      write.csv(res_pcsf()$enrich, file, row.names = F)
    }
  )
  #stops the app when browser window is closed
  session$onSessionEnded(stopApp)
}





