#' Launch Shiny App for miRNARefine
#'
#' A function that launches the Shiny app for miRNARefine.
#' The shiny app enables users to upload their own miRNA expression datasets,
#' preprocess the data according to users' selection, calculate stability
#' metrics (e.g., CV), and visualize the results. The app highlights the most
#' and least stable miRNAs on the scatter plot.
#'
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' miRNARefine::runmiRNARefine()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runmiRNARefine <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "miRNARefine")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")

  return(actionShiny)
}


# [END]
