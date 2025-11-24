#' Adaptive Filtering of miRNA Data
#'
#' Automatically filters miRNAs based on expression, variance, and missing
#' values. Thresholds are adapted to the dataset's characteristics to retain
#' informative features while removing noise.
#'
#' @param miRNAdata A data frame or matrix of miRNA expression, samples in rows, miRNAs in columns.
#' @param min_expression Minimum mean expression threshold. If NULL, use the 25th percentile of means.
#' @param min_variance Minimum variance threshold. If NULL, use the 25th percentile of variances.
#' @param max_na Maximum proportion of missing values allowed per miRNA. Default 0.2.
#' @param report_summary Logical, whether to print filtering summary. Default TRUE.
#'
#' @return A filtered data frame with informative miRNAs retained.
#'
#' @examples
#' # Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq1)
#' filtered_data <- adaptiveFiltering(miRNAdata = miRNASeq1,
#'                                    min_expression = NULL,
#'                                    min_variance = NULL,
#'                                    max_na = 0.2,
#'                                    report_summary = TRUE)
#' filtered_data
#'
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' \dontrun{
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   subset <- RTCGA.miRNASeq::ACC.miRNASeq[1:200, 1:20]
#'   sample <- subset[subset$miRNA_ID == "read_count", ]
#'   filtered <- adaptiveFiltering(miRNAdata = sample,
#'                                 min_expression = NULL,
#'                                 min_variance = NULL,
#'                                 max_na = 0.2,
#'                                 report_summary = TRUE)
#'   head(filtered)
#' }}
#'
#' @references
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Castelluzzo M, Detassis S, Denti M, Ricci L (2023).
#' \emph{MiRNAQCD: Micro-RNA Quality Control and Diagnosis}.
#' R package version 1.1.3.
#' \href{https://CRAN.R-project.org/package=MiRNAQCD}{Link}.
#'
#' Hastie, T., Tibshirani, R., & Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' @export
#' @import stats
#'
adaptiveFiltering <- function(miRNAdata,
                              min_expression = NULL,
                              min_variance = NULL,
                              max_na = 0.2,
                              report_summary = TRUE) {

  # Check if input is a matrix or dataframe
  miRNAdata <- inputCheckGeneral(miRNAdata)

  # # Check if input is a matrix or dataframe
  # if (is.matrix(miRNAdata)) {
  #   miRNAdata <- as.data.frame(miRNAdata)
  # } else if (is.data.frame(miRNAdata)) {
  #   class(miRNAdata) <- "data.frame"
  # } else {
  #   stop("`miRNAdata` must be a data frame or matrix.")
  # }
  #
  # # Stop when input is an empty dataframe
  # if (nrow(miRNAdata) == 0 || ncol(miRNAdata) == 0) {
  #   stop("Empty dataframe input")
  # }
  #
  # # Converse non-numeric data
  # miRNAdata <- as.data.frame(
  #   lapply(miRNAdata, function(x) {
  #     if (is.factor(x)) as.numeric(as.character(x)) else x
  #   })
  # )

  # Calculate missing value proportion per miRNA
  na_prop <- colSums(is.na(miRNAdata)) / nrow(miRNAdata)

  # Calculate mean and variance per miRNA
  col_mean <- colMeans(miRNAdata, na.rm = TRUE)
  col_var  <- apply(miRNAdata, 2, stats::var, na.rm = TRUE)

  # Set adaptive thresholds if not provided
  if (is.null(min_expression)) {
    min_expression <- stats::quantile(col_mean, 0.25, na.rm = TRUE)
  }
  if (is.null(min_variance)) {
    min_variance   <- stats::quantile(col_var, 0.25, na.rm = TRUE)
  }

  # Identify miRNAs to keep
  keep <- (col_mean >= min_expression) &
    (col_var >= min_variance) &
    (na_prop <= max_na)

  filtered <- miRNAdata[, keep, drop = FALSE]

  # Print summary if needed
  if (isTRUE(report_summary)) {
    message(sprintf("Original miRNAs: %d", ncol(miRNAdata)))
    message(sprintf("Filtered miRNAs: %d", ncol(filtered)))
    message(sprintf("Removed %d miRNAs", ncol(miRNAdata) - ncol(filtered)))
  }

  return(filtered)

}

