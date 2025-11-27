library(shiny)
library(shinyalert)

Sys.setenv(LANGUAGE = "en")

# Define UI ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1("miRNA Data Explorer: Preprocessing, Analysis, and
                     Stability Visualization")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Description: This is a Shiny App that is part of the miRNARefine
             R package. This Shiny application provides an interactive interface
             to explore, preprocess, and assess miRNA expression datasets. It
             integrates the core functions of the package to help researchers
             improve data quality beyond standard QC."),

      tags$ul(

        tags$li("Perform ", tags$b("adaptive filtering"),
                " to remove uninformative or noisy miRNAs."),
        tags$li("Detect and correct ", tags$b("outliers"),
                " using PCA-based diagnostics."),
        tags$li("Handle ", tags$b("missing values"),
                " with adaptive imputation methods (mean, median, KNN)."),
        tags$li("Compare and apply ", tags$b("normalization strategies"),
                " for consistent expression data."),
        tags$li("Detect and correct ", tags$b("batch effects"),
                " to avoid technical confounding."),
        tags$li("Evaluate ", tags$b("miRNA stability"),
                " using metrics like CV or MAD, identifying the most reliable
                features."),
        tags$li("Generate ", tags$b("diagnostic visualizations"),
                ", including stability distribution plots highlighting the most
                and least stable miRNAs.")
      ),

      br(),
      br(),

      # input
      tags$b("Instruction: You can either upload your own miRNA expression
             dataset in CSV format,",
             "or use the example dataset included in this package.",
             "If you upload your own CSV file,",
             "each row should correspond to a sample,",
             "and each column should correspond to a miRNA.",
             "Make sure the first row contains column headers (miRNA names)."),

      radioButtons("dataset_source", "Choose dataset:",
                   choices = c("Example dataset", "Upload your own file")),

      conditionalPanel(
        condition = "input.dataset_source == 'Upload your own file'",
        fileInput("file", "Upload CSV file", accept = ".csv")
      ),

      conditionalPanel(
        condition = "input.dataset_source == 'Example dataset'",
        radioButtons("example_choice", "Select example dataset:",
                     choices = c("Example 1", "Example 2"))
      ),

      br(),

      # Preprocessing options
      tags$h3("Preprocessing options"),
      tags$h4("0. Remove zeros"),

      # 0. Remove all-zero columns
      checkboxInput("remove_zeros", "Remove columns that are all zeros?",
                    value = FALSE),

      # action button
      actionButton("run_btn", "Run")


    ), # end of sidebar panel

    mainPanel(

      tabsetPanel(
        id = "preprocess_tabs",

        tabPanel("Raw Data Preview", tableOutput("raw_preview")),

        tabPanel("Remove all-zero columns", uiOutput("zeros_removed_preview"))
      )

    ) # end of main panel
  )
)


# Define server ----
server <- function(input, output) {

  # Reactive expression to read dataset based on user selection

  miRNAdata <- reactive({

    data <- if (input$dataset_source == "Example dataset") {

      # Load package-internal example data
      if (input$example_choice == "Example 1") {
        data("miRNASeq1", package = "miRNARefine")  # Load Example 1
        miRNASeq1
      } else {
        data("miRNASeq2", package = "miRNARefine")  # Load Example 2
        miRNASeq2
      }

    } else {
      # Read user-uploaded CSV file
      req(input$file)  # Ensure the file has been uploaded
      read.csv(input$file$datapath)
    }

    # ensure all entries are numeric
    data[] <- lapply(data, function(x) {
      if (is.factor(x)) as.numeric(as.character(x)) else x
    })

    data

  })

  # raw data preview

  output$raw_preview <- renderTable({
    head(miRNAdata())
  })


  # 0. Remove all-zero columns

  removed_zeros_data <- eventReactive(input$run_btn,{

    req(miRNAdata())
    data <- miRNAdata()

    if (isTRUE(input$remove_zeros)) {
      data <- data[, colSums(data != 0, na.rm = TRUE) > 0]
    }
    data
  })

  output$zeros_removed_preview <- renderUI({

    req(removed_zeros_data())
    if (!isTRUE(input$remove_zeros)) {
      HTML("<i>User preferred not to remove all-zero columns.</i>")

    } else {
      tableOutput("zeros_removed_table")
    }
  })

  output$zeros_removed_table <- renderTable({

    head(removed_zeros_data())

  })

}

# Create shiny app ----
shiny::shinyApp(ui = ui, server = server)


# [END]
