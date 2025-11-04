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
#' @examples
#' # Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq2)
#' filled2 <- missingValueHandling(miRNAdata = miRNASeq2,
#'                                 method = "median",
#'                                 report_summary = FALSE)
#' result <- detectBatch(miRNAdata = filled2, correct = TRUE,
#'                       report_summary = FALSE)
#' result$corrected_data
#' result$batch_effects
#'
#' \dontrun{
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'   result <- detectBatch(miRNAdata = sample, correct = TRUE,
#'                         batch = rep(1:2, each = 50)
#'                         report_summary = FALSE)
#'   result$corrected_data
#'   result$batch_effects
#'}}
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

  # When no batch factor provided, set it same for each row
  if (is.null(batch)) {
    batch <- rep("Batch1", nrow(miRNAdata))
  }

  if (length(batch) != nrow(miRNAdata)) {
    stop("Length of batch must match the number of samples.")
  }

  batch <- as.factor(batch)

  # Detect batch differences
  batch_means <- aggregate(miRNAdata, by = list(Batch = batch), FUN = mean)
  batch_vars  <- aggregate(miRNAdata, by = list(Batch = batch), FUN = var)
  summary_df <- data.frame(
    Batch = levels(batch),
    MeanAcrossGenes = rowMeans(batch_means[,-1]),
    VarAcrossGenes  = rowMeans(batch_vars[,-1])
  )

  if (isTRUE(report_summary)) {
    cat("Batch effect summary (mean ± variance across miRNAs):\n")
    print(summary_df)
  }

  # Optional correction
  corrected_data <- NULL

  if (isTRUE(correct)) {
    if (is.null(batch) || length(unique(batch)) < 2) {
      warning("Batch is NULL or has fewer than 2 levels. Skipping batch
              correction.")
      return(list(corrected_data = miRNAdata, batch_effects = summary_df))
    }

    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("Package 'sva' is required for batch correction. Please install it.")
    }
    mod <- model.matrix(~1, data = data.frame(batch))

    # Perform Batch correction
    corrected_data <- suppressMessages(sva::ComBat(dat = t(miRNAdata),
                                                  batch = batch,
                                                  mod = mod,
                                                  par.prior = TRUE))
    corrected_data <- t(corrected_data)

    if (report_summary) {
      cat("\nBatch correction applied using ComBat.\n")
    }
  }

  return(list(
    batch_effects = summary_df,
    corrected_data = corrected_data
  ))

}

