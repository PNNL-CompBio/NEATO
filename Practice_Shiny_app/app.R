library(shiny)
library(tidyverse)
library(STRINGdb)

proStat <- readRDS("/Users/lewi052/MAP/STRINGdb_exploration/proStat2.RDS")

string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")

ui <- fluidPage(selectInput(inputId = "regulation", 
                            label = "Choose whether to include upregulated or downregulated proteins in the network", 
                            choices = list("Upregulated", "Downregulated")),
                selectInput(inputId = "Comparision",
                            label = "Choose COVID-19 patient status comparision",
                            choices = list("Healthy vs Severe", "Mild vs Healthy", "Mild vs Severe")),
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
                plotOutput("network"),
                dataTableOutput("enrich")
                )


server <- function(input, output) {
  output$network <- renderPlot({
    #setting variable based on selection
    reg <- ifelse(input$regulation == "Upregulated", 1, -1)
    
    #changing names of columns based on sections
    if (input$Comparision == "Healthy vs Severe") {
      colnames(proStat)[which(names(proStat) == "P_value_A_Healthy Control_vs_Severe")] <- "P_value"
      colnames(proStat)[which(names(proStat) == "Fold_change_Healthy Control_vs_Severe")] <- "Fold_change"
      colnames(proStat)[which(names(proStat) == "Flag_A_Healthy Control_vs_Severe")] <- "Flag"
    } else if (input$Comparision == "Mild vs Healthy") {
      colnames(proStat)[which(names(proStat) == "P_value_A_Mild_vs_Healthy Control")] <- "P_value"
      colnames(proStat)[which(names(proStat) == "Fold_change_Mild_vs_Healthy Control")] <- "Fold_change"
      colnames(proStat)[which(names(proStat) == "Flag_A_Mild_vs_Healthy Control")] <- "Flag"
    } else {
      colnames(proStat)[which(names(proStat) == "P_value_A_Mild_vs_Severe")] <- "P_value"
      colnames(proStat)[which(names(proStat) == "Fold_change_Mild_vs_Severe")] <- "Fold_change"
      colnames(proStat)[which(names(proStat) == "Flag_A_Mild_vs_Severe")] <- "Flag"
    }
    
    string_db <- STRINGdb$new(version="11", species=9606, score_threshold=input$score_thresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
    
    #filtering for the appropriate regulation and significant p-value
    
    filtered_stat <- filter(proStat, Flag == reg) %>% filter(P_value < 0.05)
    
    #ordering by logFC
    filtered_stat <- filtered_stat[order(filtered_stat$`Fold_change`), ]
    
    mapped_stat <- string_db$map(filtered_stat, "Protein_Identifier", removeUnmappedRows = TRUE )
    
    hits <- mapped_stat$STRING_id[1:input$included_prots]
    wd <- getwd()
    file <- paste(wd, "hits.RDS", sep = "/")
    saveRDS(hits, file = file)
    
    string_db$plot_network(hits)
    
  })
  output$enrich <- renderDataTable({
    tmp <- input$regulation
    tmp <- input$Comparision
    tmp <- input$included_prots
    
    string_db <- STRINGdb$new(version="11", species=9606, score_threshold=input$score_thresh, input_directory="/Users/lewi052/MAP/STRINGdb_exploration/STRINGdb_cache")
    wd <- getwd()
    file <- paste(wd, "hits.RDS", sep = "/")
    hits <- readRDS(file)
    enrichment <- string_db$get_enrichment(hits)
    enrichment
  })
}

shinyApp(ui = ui, server = server)





