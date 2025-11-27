library(shiny)
library(shinyalert)

Sys.setenv(LANGUAGE = "en")
Sys.setlocale(category = "LC_ALL", locale = "C")

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

        tags$li("Handle ", tags$b("missing values"),
                " with adaptive imputation methods (mean, median, KNN)."),
        tags$li("Detect and correct ", tags$b("outliers"),
                " using PCA-based diagnostics."),
        tags$li("Perform ", tags$b("adaptive filtering"),
                " to remove uninformative or noisy miRNAs."),
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


      tags$h4("1. Missing Value Handling"),

      # 1. Missing Value Handling
      checkboxInput("do_missing_val", "Perform missing value handling?",
                    value = FALSE),

      conditionalPanel(
        condition = "input.do_missing_val == true",

        # Numeric inputs for missing value handling parameters
        radioButtons(inputId = "impute_method",
                     label = "imputation method:",
                     choices = c("median", "mean", "knn"),
                     selected = "median",
                     inline = TRUE),

        numericInput("k_value", "k value for knn imputation:",
                     value = 5, min = 1, step = 1),

        # True/False for report summary
        checkboxInput("impute_report", "Report summary?", value = TRUE)
      ),


      # action button
      actionButton("run_btn", "Run")


    ), # end of sidebar panel

    mainPanel(

      style = "max-width: 1000px; overflow-x: auto;",

      tabsetPanel(
        id = "preprocess_tabs",

        tabPanel("Raw Data Preview", tableOutput("raw_preview")),

        tabPanel("Remove all-zero columns", uiOutput("zeros_removed_preview")),

        tabPanel("Missing Value Handling", uiOutput("missing_val_preview"))
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


  # 1. Missing Value Handling

  after_imputation_result <- eventReactive(input$run_btn, {

    req(removed_zeros_data())

    tryCatch({

      msgs <- character(0)

      result <- withCallingHandlers(

        missingValueHandling(removed_zeros_data(),
                             method = input$impute_method,
                             k = input$k_value,
                             report_summary = input$impute_report),

        message = function(m) {
          msgs <<- c(msgs, m$message)
          invokeRestart("muffleMessage")
        }
      )

      list(data = result, summary = paste(msgs, collapse = "\n"))

    },error = function(e) {

      # show error message to user
      message("Error: ", e$message)

      showModal(modalDialog(
        title = "Error",
        paste("missing value handling failed:", e$message),
        easyClose = TRUE
      ))

      # stop current session
      stopApp()
    })
  })

  output$missing_val_preview <- renderUI({

    req(after_imputation_result)

    if (!isTRUE(input$do_missing_val)) {
      HTML("<i>User preferred not to perform missing value handling.</i>")

    } else {
      tagList(
        verbatimTextOutput("impute_summary"),
        tableOutput("after_imputation_table")
      )
    }
  })

  output$impute_summary <- renderText({
    req(after_imputation_result())
    after_imputation_result()$summary
  })

  output$after_imputation_table <- renderTable({
    req(after_imputation_result())
    head(after_imputation_result()$data)
  })

}

# Create shiny app ----
shiny::shinyApp(ui = ui, server = server)


# [END]
