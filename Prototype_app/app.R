library(shiny)
library(shinyjs)
library(tidyverse)
library(STRINGdb)
library(visNetwork)
library(igraph)
library(stringr)
library(shinycssloaders)
library(shinydashboard)

# defines elements of page
ui <- dashboardPage(
  dashboardHeader(title = "MAP Enrichment"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("STRING", tabName = "STRING"),
      menuItem("PCSF", tabName = "PCSF")
    )
  ),
  dashboardBody(
  tabItems(
    tabItem(tabName = "STRING",
  #separates page into sidebar and main page
  navbarPage("MAP Enrichment Analysis",
             tabPanel("Data Upload",
    #sidebar contents (inputs)
    div( id ="Sidebar",sidebarPanel(
                fileInput(inputId = "in_file",
                          label = "Upload file to be analysed",
                          accept = c(".csv",".rds"),
                          placeholder = ".csv or .rds file"),
                textInput(inputId = "species",
                          label = "Enter NCBI species identifier (e.g. Human is 9606). Identifier can be found here: https://www.ncbi.nlm.nih.gov/taxonomy",
                          value = 9606),
                textInput(inputId = "Protein",
                          label = "Enter column name containing the Protein Identifiers",
                          value = "Protein_Identifier"),
                textInput(inputId = "P_val",
                          label = "Enter column name containing P-values",
                          value = "P_value_A_Healthy Control_vs_Severe"),
                textInput(inputId = "logFC",
                          label = "Enter column name containing Log Fold Change",
                          value = "Fold_change_Healthy Control_vs_Severe"),)),
            mainPanel(
              dataTableOutput("previewTable")
            )),
    tabPanel("Network and Enrichment", 
             sidebarPanel(textInput(inputId = "included_prots",
                                      label = "Choose the number of proteins to include in the plot",
                                      value = 150),
                          textInput(inputId = "score_thresh",
                                      label = "Choose a score threshold for protein interactions",
                                      value = 400),
                          actionButton(inputId = "submit",
                                       label = "Submit Data"), width = 2),
              mainPanel(# actionButton("toggleSidebar", "Toggle sidebar"),
                        withSpinner(visNetworkOutput("network", height = "700px")),
                        withSpinner(dataTableOutput("enrich")),
                        uiOutput("downButton")
                    # tabPanel("Clusters",
                    #          downloadButton("downloadClusters", "Download Plots"),
                    #          withSpinner(plotOutput("clusters", height = "900px")))
        )
      ),
  inverse = T)),
  tabItem(tabName = "PCSF",
          actionButton(inputId = "pcsf_submit", label = "Submit"),
          withSpinner(visNetworkOutput("pcsf", height = "700px"))
          )
  ))
)

#generates output based on input
server <- function(input, output, session) {
  observeEvent(input$in_file, {
  output$previewTable <- renderDataTable({
    file <- input$in_file
    ext <- tools::file_ext(file$datapath)
    if(ext == "csv"){
      user_data <- read.csv(file$datapath)
    } else {
      user_data <- readRDS(file$datapath)
    }
    user_data
    }, options = list(pageLength = 10))})
  #program stores output of "reactive" functions so they don't have to reload when multiple elements call them
  string_db <- reactive({
    #creating stringdb object
    isolate(score_thresh <- as.numeric(input$score_thresh))
    isolate(species <- as.numeric(input$species))
    STRINGdb$new(version="11.5", species=species, score_threshold=score_thresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
  })
  hits <- reactive({
    #reading file and mapping proteins (on button press)
    if (input$submit >= 1) {
    string_db <- string_db()
    
    #isolate makes it so element isn't reloaded when values change
    isolate(file <- input$in_file)
    isolate(included_prots <- input$included_prots)
    isolate(Protein <- as.numeric(input$Protein))
    isolate(logFC <- input$logFC)
    isolate(P_val <- input$P_val)
    
    #checking file extention
    ext <- tools::file_ext(file$datapath)
    if(ext == "csv"){
      user_data <- read.csv(file$datapath)
    } else {
      user_data <- readRDS(file$datapath)
    }
    
    #converting given column names to known ones
    colnames(user_data)[which(names(user_data) == Protein)] <- "Protein_Identifier"
    colnames(user_data)[which(names(user_data) == logFC)] <- "Fold_change"
    colnames(user_data)[which(names(user_data) == P_val)] <- "P_value"
    
    #mapping based on UNIPROT identifier
    mapped_stat <- string_db$map(user_data, "Protein_Identifier", removeUnmappedRows = TRUE )
    
    #adding color based on up/down regulation
    # mapped_stat_pval05 <- string_db$add_diff_exp_color( subset(mapped_stat, P_value<0.05), logFcColStr="Fold_change" )
    #sorting by absolute value of logFC
    mapped_stat <- mutate(mapped_stat, Fold_change_abs = abs(mapped_stat$Fold_change))
    mapped_stat <- mapped_stat[order(mapped_stat$`Fold_change_abs`, decreasing = T), ]
    
    nodes <- mapped_stat[1: included_prots, c("Protein_Identifier", "STRING_id", "Gene_Name")]
    colnames(nodes)  <- c("label", "id", "name")
    edges <- string_db$get_interactions(mapped_stat$STRING_id)
    edges <- edges[!duplicated(edges),]
    colnames(edges) <- c("from", "to", "width")
    edges$width <- as.integer(as.numeric(edges$width)/100)
    print(head(edges))
    list(nodes = nodes,
         edges = edges[,1:2])
  }
  })

#NETWORK PLOTS
  #clustering nodes
  nodes <- reactive({
    string_db <- string_db()
    hits <- hits()
    nodes <- hits$nodes
    
    clusterList <- string_db$get_clusters(nodes$id)
    
    clusters <- c()
    for (x in nodes$id) {
      clusters <- c(clusters, grep(x, clusterList))
    }
    nodes$group <- clusters
    
    nodes <- nodes[!duplicated(nodes$id),]
    nodes <- nodes[order(nodes$group), ]
    for (y in nodes$group){
      if (sum(nodes$group==y) == 1){
        nodes$group[which(nodes$group == y)] <- "No Cluster"
      }
    }
    nodes
  })
  output$network <- renderVisNetwork({
    #creating network (on button press)
    if(input$submit >= 1) {
      nodes <- nodes()
      colnames(nodes) <- c("x", "id", "label", "group")
      # nodes <- nodes[!duplicated(nodes$label), ]
      hits <- hits()
      edges <- hits$edges
      # edges <- edges[nodes$id %in% edges$to & nodes$id %in% edges$from, ]
      # graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
      
      # visIgraph(graph, physics = F)
      visNetwork(nodes = nodes, edges = edges) %>% 
        visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
        visEdges(width = 10) %>%
        visInteraction(navigationButtons = T) %>% visIgraphLayout(randomSeed = 123) %>% visLegend(zoom = F) %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('sel_node', nodes.nodes);
                ;}") %>%
        visExport(type = "png", name = "exported-network", float = "right",
                  label = "Export PNG", background = "white", style= "")
    }
  })

  observeEvent(input$sel_node,{
    nodes <- nodes()
    sel_node_row <- nodes[nodes$id == input$sel_node, ]
    node_cluster <- filter(nodes, group == sel_node_row$group)
    # node_exclude <- filter(nodes, group != sel_node_row$group)
    visNetworkProxy("network") %>% visSelectNodes(node_cluster$id) %>% visFocus(id = input$sel_node, scale = 0.4)
  })
  
#ENRICHMENT TABLE
  enrich_table <- reactive({
    string_db <- string_db()
    nodes <- nodes()
    print(input$sel_node)
    #filters enrichment table based on selected cluster
    if (!is.null(input$sel_node)) {
      sel_node_row <- nodes[nodes$id == input$sel_node, ]
      nodes <- filter(nodes, group == sel_node_row$group)
    }
    enrichment <- string_db$get_enrichment(nodes$id)
    enrichment
  })
  output$enrich <- renderDataTable({
    if(input$submit >= 1) {
    #creates enrichment table
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
  output$pcsf <- renderVisNetwork({
    if (input$pcsf_submit >= 1) {
    
    isolate(file <- input$in_file)
    ext <- tools::file_ext(file$datapath)
    if(ext == "csv"){
      user_data <- read.csv(file$datapath)
    } else {
      user_data <- readRDS(file$datapath)
    }
    proStat_HvS <- user_data %>% filter(`Flag_A_Healthy Control_vs_Severe` != 0)
    
    geneList_HvS <- abs(proStat_HvS[,10])
    gene_sym <- UniProt.ws::select(org.Hs.eg.db, proStat_HvS[,17], "SYMBOL","UNIPROT")
    gene_sym <- gene_sym[!duplicated(gene_sym[,1]),]
    names(geneList_HvS) = as.character(gene_sym[,2])
    geneList_HvS <- sort(geneList_HvS, decreasing = T)
    terminals <- geneList_HvS[1:50]
    ppi <- construct_interactome(STRING)
    subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
    visIgraph(subnet) %>% visNodes(value = 45, font = list(size = 40), scaling = list(max = 75), shadow = list(enabled = T, size = 10)) %>%
      visInteraction(navigationButtons = T) %>%
      visExport(type = "png", name = "exported-network", float = "right",
                                    label = "Export PNG", background = "white", style= "")
  }})
  
  #stops the app when browser window is closed
  session$onSessionEnded(stopApp)
}

shinyApp(ui = ui, server = server)





