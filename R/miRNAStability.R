#' miRNA Stability Analysis
#'
#' Calculates feature-level stability metrics such as Coefficient of Variation
#'  (CV) and Median Absolute Deviation (MAD) for each miRNA across samples.
#' Helps identify the most stable miRNAs for reliable downstream analyses.
#'
#'
#' @param miRNAdata A numeric matrix or data frame with miRNAs as columns and samples as rows.
#' @param metrics A character vector specifying which stability metrics to calculate.
#'        Options: \code{c("CV", "MAD")}. Default is both.
#' @param report_summary Logical, whether to print a short summary of the most/least stable miRNAs. Default TRUE.
#'
#' @return A list with:
#' \describe{
#'   \item{stability_scores}{A data frame containing CV and/or MAD for each miRNA.}
#'   \item{most_stable}{Names of the top 5 most stable miRNAs.}
#'   \item{least_stable}{Names of the top 5 least stable miRNAs.}
#' }
#'
#' @examples
#' # Example 1:
#' # Using miRNASeq1 available with package
#' data(miRNASeq1)
#' result <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)
#'
#' result$most_stable
#' result$least_stable
#'
#' #Example 2:
#' # Using miRNASeq2 available with package
#' data(miRNASeq2)
#' result2 <- miRNAStability(miRNAdata = miRNASeq2, metrics = "CV",
#'                           report_summary = FALSE)
#' result2$most_stable
#' result2$least_stable
#'
#' \dontrun{
#' # Example 3:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'   result <- miRNAStability(miRNAdata = sample, report_summary = FALSE)
#'
#'   result3$most_stable
#'   result3$least_stable
#'}}
#'
#' @references
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Hastie, T., Tibshirani, R., & Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' Hima Bindu, A., Suresh, P., Reddy, R., et al. (2019). Systematic evaluation
#' of miRNA stability across biological replicates. \emph{Frontiers in Molecular
#' Biosciences}, 6, 153.
#'
#' @export
#' @import tidyr
#' @importFrom dplyr group_by summarise arrange
#'
miRNAStability <- function(miRNAdata,
                           metrics = c("CV", "MAD"),
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

  metrics <- match.arg(metrics, choices = c("CV", "MAD"), several.ok = TRUE)

  # Reshape and compute metrics using dplyr
  df_long <- miRNAdata %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "miRNA",
                        values_to = "Expression")

  results <- df_long %>%
    dplyr::group_by(miRNA) %>%
    dplyr::summarise(
      CV = if ("CV" %in% metrics) sd(Expression) / mean(Expression)
           else NA_real_,
      MAD = if ("MAD" %in% metrics) mad(Expression)
            else NA_real_,
      .groups = "drop"
    )

  # Identify stable/unstable miRNAs
  stability_metric <- if ("CV" %in% metrics) "CV" else "MAD"

  stability_order <- results %>%
    dplyr::arrange(.data[[stability_metric]])

  most_stable <- head(stability_order$miRNA, 5)
  least_stable <- tail(stability_order$miRNA, 5)

  # Optional summary
  if (isTRUE(report_summary)) {
    cat("Most stable miRNAs (lowest variability):\n")
    print(most_stable)
    cat("\nLeast stable miRNAs (highest variability):\n")
    print(least_stable)
  }

  return(list(
    stability_scores = results,
    most_stable = most_stable,
    least_stable = least_stable
  ))

}

