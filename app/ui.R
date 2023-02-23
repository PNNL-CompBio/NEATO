# defines elements of page
ui <- navbarPage("NEATO",

                 #PAGE 0: About
                 tabPanel("About",
                          column(3),
                          column(6,
                                 shiny::includeMarkdown("www/about_neato.md")
                          ),
                          column(3)
                          # column(6,
                          #   img(src = "www/neato.jpeg",
                          #       height = 1000,
                          #       width = 1000)
                          # )
                 ),

                 # PAGE 1: Data Upload
                 tabPanel("Data Upload",
                          column(3,
                                 bsCollapse(
                                   open = "Data Upload",
                                   id = "data_upload",
                                   bsCollapsePanel("Data Upload",
                                                   fluidRow(
                                                   column(9,
                                                   fileInput(inputId = "inFile",
                                                             label = "Upload Proteomics Data",
                                                             accept = c(".csv", ".tsv", ".rds"),
                                                             placeholder = ".csv, .tsv or .rds file")),
                                                   column(3,
                                                   actionButton("fileHelp", "", icon = icon("question"))
                                                   )),
                                                   fluidRow(
                                                   column(9,
                                                   checkboxInput(inputId = "testDataButton", label = "Load Test Data"),
                                                   ),
                                                   column(3,
                                                   actionButton("testDataHelp", "", icon = icon("question"))
                                                   ))),
                                   bsCollapsePanel("Data Specifics",
                                                   selectInput(inputId = "Protein",
                                                               label = "Select column name containing Protein Identifiers",
                                                               # value = "Reference",
                                                               choices = c("")
                                                   ),
                                                   selectInput(inputId = "logFC",
                                                               label = "Select column name containing Log-Fold Changes",
                                                               # value = "Fold_change_Mock_vs_Infection"),
                                                               choices = c("")
                                                   ),
                                                   selectInput(inputId = "P_val",
                                                               label = "Select column name containing P-values",
                                                               # value = "P_value_A_Mock_vs_Infection"
                                                               choices = c("")
                                                   ),
                                                   selectInput(inputId = "species",
                                                               label = "Select the species you wish to analyze",
                                                               choices = c("Saccharomyces cerevisiae (Yeast)" = 4932,
                                                                           "Mus musculus (Mouse)" = 10090,
                                                                           "Homo Sapien (Human)" = 9606,
                                                                           "Arabidopsis thaliana (Thale cress)" = 3702))),
                                   bsCollapsePanel("Parameters",

                                                   selectInput(inputId = "algorithm",
                                                               label = "Select algorithm to generate graph",
                                                               choices = c("Explicit graph (only use proteins in uploaded data)" = "explicitGraph",
                                                                           "Inferred graph (include proteins not in uploaded data" = "inferredGraph")),
                                                   textInput(inputId = "P_val_cut",
                                                             label = "Set P-value cutoff for your data",
                                                             value = 0.05),
                                                   textInput(inputId = "includedProts",
                                                             label = "Choose the number of proteins to include in the analysis",
                                                             value = 150),
                                                   textInput(inputId = "scoreThresh",
                                                             label = "Choose a score threshold for protein interactions",
                                                             value = 500)
                                   ),
                                   bsCollapsePanel("Sumbit Job",
                                                   fluidRow(
                                                     actionButton(inputId = "sendService", label = "Process Uploaded Data"),
                                                     actionButton(inputId = "checkStatus", label = "Check Status")),
                                                     # withSpinner(textOutput("Loading")),
                                                     textOutput("Submitted")
                                   ))),

                          column(9,
                                 bsCollapse(
                                   id = "data_preview",
                                   open = "Data Preview",
                                   bsCollapsePanel("Data Preview",
                                                   dataTableOutput("previewTable")
                                   )))),

                 # PAGE 2: OUTPUT
                 tabPanel("Output",
                          actionButton(inputId = "pcsfSubmit", label = "Display Output"),
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
                 # icon = "www/neato.png",
                 # theme = "bootstrap.css",
                 selected = "Data Upload",
                 inverse = T)