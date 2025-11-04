#' Compare Normalization Methods for miRNA Data
#'
#' Applies multiple normalization methods (log2 transformation, z-score scaling,
#' and quantile normalization) and compares their effects on miRNA expression
#' data. Helps select the most appropriate normalization strategy for consistent
#' and comparable expression levels.
#'
#' @param miRNAdata A numeric data frame or matrix, with samples in rows and miRNAs in columns.
#' @param methods A character vector of normalization methods to apply.
#' Supported: "log2", "zscore", "quantile". Default is all three.
#' @param report_summary Logical, whether to print summary information. Default TRUE.
#' @param choose_best Logical, whether to automatically recommend the best method. Default FALSE.
#'
#' @return A result list containing:
#'   \itemize{
#'     \item normalized: list of normalized datasets
#'     \item best_method: recommended method if choose_best = TRUE
#'   }
#'
#' @examples
#' Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq1)
#' result <- compareNormalization(miRNAdata = miRNASeq1,
#'                                choose_best = TRUE)
#' result$normalized$log2[1:5, 1:5]
#' result$best_method
#'
#' #' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' \dontrun{
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'   result <- compareNormalization(miRNAdata = sample,
#'                                  methods = c("zscore", "quantile"),
#'                                  report_summary = TRUE,
#'                                  choose_best = TRUE)
#'   head(filtered)
#' }}
#'
#' @references
#' Bolstad, B. M., Irizarry, R. A., Astrand, M., & Speed, T. P. (2003).
#' A comparison of normalization methods for high density oligonucleotide
#' array data based on variance and bias. Bioinformatics, 19(2), 185–193.
#'
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Castelluzzo M, Detassis S, Denti M, Ricci L (2023).
#' \emph{MiRNAQCD: Micro-RNA Quality Control and Diagnosis}.
#' R package version 1.1.3.
#' \href{https://CRAN.R-project.org/package=MiRNAQCD}{Link}.
#'
#' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Hastie, T., Tibshirani, R., & Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' @export
#' @import preprocessCore
#' @importFrom stats var
#'
compareNormalization <- function(miRNAdata,
                                 methods = c("log2", "zscore", "quantile"),
                                 report_summary = TRUE,
                                 choose_best = FALSE) {

  # Check if input is valid
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

  # Match methods
  methods <- match.arg(methods, several.ok = TRUE)

  # Initialize result list
  normalized <- list()

  # Log2 transformation
  if ("log2" %in% methods) {
    # Add small constant to avoid log2(0)
    normalized$log2 <- log2(miRNAdata + 1)
  }

  # Z-score scaling (per column)
  if ("zscore" %in% methods) {
    normalized$zscore <- as.data.frame(stats::scale(miRNAdata))
  }

  # Quantile normalization
  if ("quantile" %in% methods) {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("The 'preprocessCore' package is required for quantile
           normalization.")
    }
    normalized$quantile <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(miRNAdata)))
    colnames(normalized$quantile) <- colnames(miRNAdata)
    rownames(normalized$quantile) <- rownames(miRNAdata)
  }

  # Optional summary
  if (isTRUE(report_summary)) {
    message("Normalization applied: ", paste(names(normalized), collapse = ", "))
    for (nm in names(normalized)) {
      message(sprintf("%s: min=%g, max=%g", nm,
                      min(normalized[[nm]], na.rm = TRUE),
                      max(normalized[[nm]], na.rm = TRUE)))
    }
  }

  # Choose best method based on variance uniformity
  best_method <- NULL
  if (isTRUE(choose_best)) {
    score <- sapply(normalized, function(x) {
      vars <- apply(x, 2, stats::var, na.rm = TRUE)
      sd(vars)  # smaller SD of variances = more uniform
    })
    best_method <- names(which.min(score))
    if (isTRUE(report_summary)) {
      message(sprintf("Recommended method: %s", best_method))
    }
  }

  return(list(normalized = normalized, best_method = best_method))

}
