# defines elements of page
ui <- navbarPage("NEATO",

      # PAGE 1: Data Upload
      tabPanel("Data_Upload",
              column(6,
              bsCollapse(
                open = "Data Upload",
                id = "data_upload",
                bsCollapsePanel("Data Upload",
                  column(6,
                  fileInput(inputId = "inFile",
                            label = "Upload file to be analysed",
                            accept = c(".csv",".rds"),
                            placeholder = ".csv or .rds file"),
                  textInput(inputId = "Protein",
                            label = "Enter column name containing the Protein Identifiers",
                            value = "Protein_Identifier"
                            ),
                  textInput(inputId = "logFC",
                            label = "Enter column name containing Log Fold Change",
                            value = "Fold_change_Healthy Control_vs_Severe"
                            ),
                  textInput(inputId = "P_val",
                            label = "Enter column name containing P-values",
                            value = "P_value_A_Healthy Control_vs_Severe"
                            ),
                  fluidRow(
                  actionButton(inputId = "sendService", label = "Process Uploaded Data"),
                  actionButton(inputId = "checkStatus", label = "Check Status"))),
                  column(6,
                  textInput(inputId = "species",
                           label = "Enter NCBI species identifier (e.g. Human is 9606).",
                           value = 9606),
                  selectInput(inputId = "algorithm",
                              label = "Choose algorithm to generate graph",
                              choices = c("Explicit graph (only use proteins in uploaded data)" = "explicitGraph",
                                          "Inferred graph (include proteins not in uploaded data" = "inferredGraph")),
                  textInput(inputId = "P_val_cut",
                            label = "Choose P-value cutoff for the data",
                            value = 0.05),
                  textInput(inputId = "includedProts",
                            label = "Choose the number of proteins to include in the plot",
                            value = 150),
                  textInput(inputId = "scoreThresh",
                            label = "Choose a score threshold for protein interactions",
                            value = 500),
                  )))),

              column(6,
                bsCollapse(
                  id = "data_preview",
                  open = "Data Preview",
                  bsCollapsePanel("Data Preview",
                    dataTableOutput("previewTable")
                                    )))),

      # PAGE 2: OUTPUT
      tabPanel("Output",
              actionButton(inputId = "pcsfSubmit", label = "Generate Plot"),
              bsCollapse(
                open = c("Network", "Enrichment Table"),
                multiple = T,
                id = "pcsf_net",
                bsCollapsePanel("Network",
                  withSpinner(visNetworkOutput("pcsf", height = "700px"))),
                bsCollapsePanel("Enrichment Table",
                  splitLayout(cellWidths = c("30%", "30%", "40%"),
                  textInput(inputId = "pvalFilterPCSF", label = "Enter a P-Value Cutoff", value = 1),
                  textInput(inputId = "BHpvalFilterPCSF", label = "Enter a BH P-Value Cutoff", value = 1),
                  radioGroupButtons(inputId = "bgChoicePCSF", 
                                    label = "Use Full Dataset of Network as the Background for Enrichent",
                                    choices = c("Full Dataset" = "full", "Network" = "net"),
                                    selected = "full"
                  )),
                  withSpinner(dataTableOutput("pcsf_enrich")),
                  uiOutput("downPCSFButtonFull"),
                  uiOutput("downPCSFButtonCurrent")))
      ),
    # # PAGE 4: SPRAS
    # tabItem(tabName = "SPRAS_test",
    #         actionButton(inputId = "sprasSubmit", label = "Generate Files"),
    #         withSpinner(visNetworkOutput("sprasNet", height = "700px"))
    # )

# theme = "bootstrap.css"
inverse = T)
