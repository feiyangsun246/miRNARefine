library(shiny)
library(shinyalert)

# Define UI ----
ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1("miRNA Data Explorer: Preprocessing, Analysis, and Stability Visualization")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Description: This is a Shiny App that is part of the miRNARefine
             R package. "),

      br()

    ), # end of Sidebar panel

    mainPanel(

    )
  )
)

# Define server ----
server <- function(input, output) {
}

# Create shiny app ----
shiny::shinyApp(ui = ui, server = server)


# [END]
