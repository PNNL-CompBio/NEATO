library(shiny)
library(tidyverse)
library(STRINGdb)
library(shinycssloaders)

# defines elements of page
ui <- fluidPage(
  #separates page into sidebar and main page
  sidebarLayout(
    #sidebar contents (inputs)
    sidebarPanel(
                fileInput(inputId = "in_file",
                          label = "Upload file to be analysed",
                          accept = c(".csv",".rds"),
                          placeholder = ".csv or .rds file"),
                textInput(inputId = "Protein",
                          label = "Enter column name containing the Protein Identifiers",
                          value = "Protein_Identifier"),
                textInput(inputId = "P_val",
                          label = "Enter column name containing p_values",
                          value = "P_value_A_Healthy Control_vs_Severe"),
                textInput(inputId = "logFC",
                          label = "Enter column name containing Log Fold Change",
                          value = "Fold_change_Healthy Control_vs_Severe"),
                sliderInput(inputId = "included_prots",
                            label = "Choose the number of proteins to include in the plot",
                            min = 10,
                            max = 500,
                            value = 150),
                sliderInput(inputId = "score_thresh",
                            label = "Choose a confidence threshold for protein interactions",
                            min = 100,
                            max = 800,
                            value = 400),
                actionButton(inputId = "submit",
                             label = "Sumbit Data"),
  ),
      #main page (outputs on tabs)
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Network Plot", 
                             withSpinner(plotOutput("network", height = "700px")),
                             uiOutput("legend")),
                    # dataTableOutput("top"),
                    tabPanel("Enrichment", withSpinner(dataTableOutput("enrich"))),
                    # plotOutput("clusters")
        )
      )
  )
)

#generates output based on input
server <- function(input, output, session) {
  #program stores output of "reactive" functions so they don't have to reload when multiple elements call them
  string_db <- reactive({
    #creating stringdb object
    isolate(score_thresh <- input$score_thresh)
    STRINGdb$new(version="11", species=9606, score_threshold=score_thresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
  })
  hits <- reactive({
    #reading file and mapping proteins (on button press)
    if(input$submit >= 1) {
    string_db <- string_db()
    
    #isolate makes it so element isn't reloaded when values change
    isolate(file <- input$in_file)
    isolate(included_prots <- input$included_prots)
    isolate(Protein <- input$Protein)
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
    mapped_stat_pval05 <- string_db$add_diff_exp_color( subset(mapped_stat, P_value<0.05), logFcColStr="Fold_change" )
    #sorting by absolute value of logFC
    mapped_stat_pval05 <- mutate(mapped_stat_pval05, Fold_change_abs = abs(mapped_stat_pval05$Fold_change))
    mapped_stat_pval05 <- mapped_stat_pval05[order(mapped_stat_pval05$`Fold_change_abs`, decreasing = T), ]
    
    #grabbing proteins with the highest absolute logFC
    mapped_stat_pval05[1:included_prots, ]
  }
  })
  output$network <- renderPlot({
    #creating network (on button press)
    if(input$submit >= 1) {
      string_db <- string_db()
      hits <- hits()
      payload_id <- string_db$post_payload( hits$STRING_id, colors=hits$color )
      string_db$plot_network(hits$STRING_id, payload_id = payload_id)
    }
  })
  output$legend <- renderUI({
    #adding legend for STRINGdb network graphs
    if(input$submit >= 1) {
      tags$iframe(style="height:610px; width:100%", src="legend.png")
    }
  })
  # output$top <- renderDataTable({
  #   #calling all input variables so this output always updates appropriately
  #   tmp <- input$Comparision
  #   tmp <- input$included_prots
  #   tmp <- input$score_thresh
  #   #reading in file
  #   wd <- getwd()
  #   file <- paste(wd, "hits.RDS", sep = "/")
  #   hits <- readRDS(file)
  #   #Posting only some columns for clarity
  #   hits[, c("Protein_Identifier", "STRING_id", "Fold_change", "Fold_change_abs", "P_value")]
  # })
  output$enrich <- renderDataTable({
    #creates enrichment table
    if(input$submit >= 1) {
    string_db <- string_db()
    hits <- hits()
    
    # string_db$set_background(user_data$Protein_Identifier)
    enrichment <- string_db$get_enrichment(hits$STRING_id)
    #cutting a couple of rows to make table fit on the screen
    enrichment[, c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr")]
    }
    })
  # output$clusters <-  renderPlot({
  #   tmp <- input$Comparision
  #   tmp <- input$included_prots
  # 
  #   string_db <- STRINGdb$new(version="11", species=9606, score_threshold=input$score_thresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
  #   wd <- getwd()
  #   file <- paste(wd, "hits.RDS", sep = "/")
  #   hits <- readRDS(file)
  #   clustersList <- string_db$get_clusters(hits$STRING_id)
  #   par(mfrow=c(2,2))
  #   for(i in seq(1:4)){
  #     string_db$plot_network(clustersList[[i]])
  #   }
  # })
  
  #stops the app when browser window is closed
  session$onSessionEnded(stopApp)
}

shinyApp(ui = ui, server = server)





