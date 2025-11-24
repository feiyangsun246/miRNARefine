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
#' \describe{
#'    \item{normalized}{list of normalized datasets}
#'    \item{best_method}{recommended method if \code{choose_best = TRUE}}
#' }
#'
#' @examples
#' # Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq1)
#' result <- compareNormalization(miRNAdata = miRNASeq1,
#'                                choose_best = TRUE)
#' result$normalized$log2[1:5, 1:5]
#' result$best_method
#'
#' \dontrun{
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset (Chodor, 2025)
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   subset <- RTCGA.miRNASeq::ACC.miRNASeq[1:200, 1:20]
#'   sample <- subset[subset$miRNA_ID == "read_count", ]
#'   result <- compareNormalization(miRNAdata = sample,
#'                                  methods = c("zscore", "quantile"),
#'                                  report_summary = TRUE,
#'                                  choose_best = TRUE)
#'   head(result)
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
#' @import preprocessCore stats
#'
compareNormalization <- function(miRNAdata,
                                 methods = c("log2", "zscore", "quantile"),
                                 report_summary = TRUE,
                                 choose_best = TRUE) {

  # Check if input is a matrix or dataframe
  miRNAdata <- inputCheckGeneral(miRNAdata)

  # Stop when there is NA in the dataset
  if (any(is.na(miRNAdata))) {
    stop("Dataset contains missing values (NA).
         Consider running missingValueHandling() first.")
  }

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
    # There might be constant columns
    constant_cols <- apply(miRNAdata, 2, sd, na.rm = TRUE) == 0
    zscore_norm <- as.data.frame(scale(miRNAdata))
    if (any(constant_cols)) {
      warning("Some columns have zero variance; z-score normalization skipped
              for these columns.")
      zscore_norm[, constant_cols] <- miRNAdata[, constant_cols]
    }
    normalized$zscore <- zscore_norm
  }

  # Quantile normalization using preprocessCore (Bolstad, 2025)
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

  # Choose best method based on variance uniformity (Castelluzzo et al., 2023; Hastie et al., 2009)
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

