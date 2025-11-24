#' Check and standardize input miRNA data
#'
#' This function validates and standardizes the input miRNA expression data.
#' It ensures the data is a matrix or data frame, is non-empty, and converts
#' factor columns to numeric.
#'
#' @param miRNAdata Input dataset (mostly a matrix or data frame) containing miRNA expression values.
#'
#' @return A data frame with numeric columns after necessary conversions.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item Ensures the input is a matrix or data frame.
#'   \item Stops with an error if the input is empty (no rows or columns).
#'   \item Converts factor columns to numeric via character conversion.
#' }
#'
#' @examples
#' # Example 1:
#' # Input a dataframe
#' df <- data.frame(a = factor(c("1", "2")), b = c(3, 4))
#' inputCheckGeneral(df)
#'
#' # Example 2:
#' # Input a matrix
#' mat <- matrix(1:4, nrow = 2)
#' inputCheckGeneral(mat)
#'
#' @export

inputCheckGeneral <- function(miRNAdata) {
  # Check if input is a matrix or dataframe
  if (is.matrix(miRNAdata)) {
    miRNAdata <- as.data.frame(miRNAdata)
  } else if (is.data.frame(miRNAdata)) {
    class(miRNAdata) <- "data.frame"
  } else {
    stop("`miRNAdata` must be a data frame or matrix.")
  }

  # Stop when input is an empty dataframe
  if (nrow(miRNAdata) == 0 || ncol(miRNAdata) == 0) {
    stop("Empty dataframe input")
  }

  # Converse non-numeric data
  miRNAdata <- as.data.frame(
    lapply(miRNAdata, function(x) {
      if (is.factor(x)) as.numeric(as.character(x)) else x
    })
  )

  return(miRNAdata)
}
