#' Detect and Correct Batch Effects in miRNA Data
#'
#' Detects and optionally corrects batch effects in the dataset. Ensures that
#' technical or experimental batch differences do not confound downstream
#' analyses.
#'
#' @param miRNAdata A numeric data frame or matrix where rows are samples and columns are miRNAs.
#' @param batch A factor or character vector indicating the batch assignment for each sample. Default NULL.
#' @param correct Logical, whether to perform batch effect correction using ComBat. Default FALSE.
#' @param report_summary Logical, whether to print summary information of batch differences before (and after) correction. Default TRUE.
#'
#' @return A list containing:
#' \describe{
#'    \item{batch_effects}{Data frame summarizing per-batch means and variances.}
#'    \item{corrected_data}{Batch-corrected expression matrix if \code{correct = TRUE}}
#' }
#'
#' @references
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' #' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' #' Hastie, T., Tibshirani, R., & Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' \emph{Biostatistics}, 8(1), 118–127.
#'
#' Leek, J. T., et al. (2010). Tackling the widespread and critical impact
#' of batch effects in high-throughput data. \emph{Nature Reviews Genetics},
#' 11(10), 733–739.
#'
#' @export
#' @importFrom stats prcomp aov
#' @import sva
#'
detectBatch <- function(miRNAdata, batch = NULL, correct = FALSE,
                        report_summary = TRUE) {

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

  # Stop when there is NA in the dataset
  if (any(is.na(miRNAdata))) {
    stop("Dataset contains missing values (NA).
         Consider running missingValueHandling() first.")
  }
}
