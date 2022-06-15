# defines elements of page
ui <- dashboardPage(
  
  # HEADER
  dashboardHeader(title = "MAP Enrichment"),
  
  # SIDEBAR
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "Data_Upload"),
      menuItem("STRING", tabName = "STRING"),
      menuItem("PCSF", tabName = "PCSF")
    )
  ),
  
  # BODY CONTENT
  dashboardBody(
    tabItems(
      
      # PAGE 1: Data Upload
      tabItem(tabName = "Data_Upload",
              box(
              fileInput(inputId = "inFile",
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
                        value = "Fold_change_Healthy Control_vs_Severe"),
              textInput(inputId = "includedProts",
                        label = "Choose the number of proteins to include in the plot",
                        value = 150),
              textInput(inputId = "scoreThresh",
                        label = "Choose a score threshold for protein interactions",
                        value = 400),
              actionButton(inputId = "construct",
                           label = "Construct Databases"),
              withSpinner(textOutput("constructed"))),
              box(
              dataTableOutput("previewTable")
                                  )),
      # PAGE 2: STRING
      tabItem(tabName = "STRING",
                  actionButton(inputId = "submit", label = "Generate Plot"),
                  withSpinner(visNetworkOutput("network", height = "700px")),
                  withSpinner(dataTableOutput("enrich")),
                  uiOutput("downButton")
    ),
      
      # PAGE 3: PCSF
      tabItem(tabName = "PCSF",
              actionButton(inputId = "pcsfSubmit", label = "Generate Plot"),
              withSpinner(visNetworkOutput("pcsf", height = "700px")),
                withSpinner(dataTableOutput("pcsf_enrich")),
              uiOutput("downPCSFButton")
      )
)))
