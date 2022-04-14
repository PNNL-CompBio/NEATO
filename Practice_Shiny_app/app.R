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
                plotOutput("network"))


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
    
    #filtering for the appropriate regulation and significant p-value
    
    filtered_stat <- filter(proStat, Flag == reg) %>% filter(P_value < 0.05)
    
    #ordering by logFC
    filtered_stat[order(filtered_stat$`Fold_change`), ]
    
    mapped_stat <- string_db$map(filtered_stat, "Protein_Identifier", removeUnmappedRows = TRUE )
    if (nrow(mapped_stat) >= 150) {
    hits <- mapped_stat$STRING_id[1:150]
    string_db$plot_network(hits)
    } else {
    hits <- mapped_stat$STRING_id[1:150]
    string_db$plot_network(hits)
    }
    
  })
}

shinyApp(ui = ui, server = server)