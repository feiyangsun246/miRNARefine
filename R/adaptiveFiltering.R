#' Adaptive Filtering of miRNA Data
#'
#' Automatically filters miRNAs based on expression, variance, and missing
#' values.
#' Thresholds are adapted to the dataset's characteristics to retain informative
#' features while removing noise.
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
#' filtered_data <- adaptiveFiltering(miRNASeq1,
#'                                    min_expression = NULL,
#'                                    min_variance = NULL,
#'                                    max_na = 0.2,
#'                                    report_summary = TRUE)
#' filtered_data
#'
#' @export
#' @import stats
#'
adaptiveFiltering <- function(miRNAdata,
                              min_expression = NULL,
                              min_variance = NULL,
                              max_na = 0.2,
                              report_summary = TRUE) {

  # Input check
  if (!is.data.frame(miRNAdata) && !is.matrix(miRNAdata)) {
    stop("`data` must be a data frame or matrix with samples in rows and miRNAs in columns.")
  }

  # Calculate missing value proportion per miRNA
  na_prop <- colSums(is.na(miRNAdata)) / nrow(miRNAdata)

  # Calculate mean and variance per miRNA
  col_mean <- colMeans(miRNAdata, na.rm = TRUE)
  col_var  <- apply(miRNAdata, 2, stats::var, na.rm = TRUE)

  # Set adaptive thresholds if not provided
  if (is.null(min_expression)) min_expression <- stats::quantile(col_mean, 0.25)
  if (is.null(min_variance))   min_variance   <- stats::quantile(col_var, 0.25)

  # Identify miRNAs to keep
  keep <- (col_mean >= min_expression) &
    (col_var >= min_variance) &
    (na_prop <= max_na)

  filtered <- miRNAdata[, keep, drop = FALSE]

  # Print summary if needed
  if (report_summary) {
    message(sprintf("Original miRNAs: %d", ncol(miRNAdata)))
    message(sprintf("Filtered miRNAs: %d", ncol(filtered)))
    message(sprintf("Removed %d miRNAs", ncol(miRNAdata) - ncol(filtered)))
  }

  filtered <- miRNAdata
  return(filtered)

}
