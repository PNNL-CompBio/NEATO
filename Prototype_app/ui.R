# defines elements of page
ui <- dashboardPage(skin = "black",

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
              column(3,
              bsCollapse(
                open = "Data Upload",
                id = "data_upload",
                bsCollapsePanel("Data Upload",
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
                            value = 500)))),
              # actionButton(inputId = "construct",
              #              label = "Construct Databases"),
              # withSpinner(textOutput("constructed"))),
              column(9,
                bsCollapse(
                  id = "data_preview",
                  open = "Data Preview",

                  bsCollapsePanel("Data Preview",
                    dataTableOutput("previewTable")
                                    )))),
      # PAGE 2: STRING
      tabItem(tabName = "STRING",
                  actionButton(inputId = "submit", label = "Generate Plot"),
              bsCollapse(
                open = c("Network", "Enrichment Table"),
                multiple = T,
                id = "net_enrich",
                bsCollapsePanel("Network",
                  withSpinner(visNetworkOutput("network", height = "700px"))),
                bsCollapsePanel("Enrichment Table",
                  splitLayout(cellWidths = c("30%", "30%", "40%"),
                  textInput(inputId = "pvalFilter", label = "Enter a P-Value Cutoff", value = 1),
                  textInput(inputId = "BHpvalFilter", label = "Enter a BH P-Value Cutoff", value = 1),
                  radioGroupButtons(inputId = "bgChoice", 
                                    label = "Use Full Dataset of Network as the Background for Enrichent",
                                    choices = c("Full Dataset" = "full", "Network" = "net"),
                                    selected = "full"
                  )),
                  withSpinner(dataTableOutput("enrich")),
                  uiOutput("downButtonFull"),
                  uiOutput("downButtonCurrent"))
    )),

      # PAGE 3: PCSF
      tabItem(tabName = "PCSF",
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
      )
)))
