#' Detect Outliers in miRNA Data Using PCA
#'
#' Identifies anomalous samples using principal component analysis (PCA)
#' following the approach of Jolliffe (2002). Samples that deviate strongly from
#' the multivariate center are flagged as outliers to improve data quality.
#'
#' @param miRNAdata A numeric matrix or data frame, samples in rows, miRNAs in columns.
#' @param n_components Number of principal components to use. Default 2.
#' @param row_threshold Quantile threshold for flagging outliers based on PCA distance. Default 0.99.
#' @param z_threshold Number of standard deviations for column-level outliers. Default 3.
#' @param whether_scale Logical, whether to standardize variables before PCA. Default TRUE.
#' @param report_summary Logical, whether to print the number of outliers detected. Default TRUE.
#'
#' @return A list with:
#'   - sample_outliers: logical vector for each sample
#'   - miRNA_outliers: logical matrix for each sample x miRNA
#'   - pca_res: prcomp object
#'   - scores: PCA scores
#'
#' @examples
#' \dontrun{
#' # Example 1:
#' Using miRNASeq data available with package
#' result <- detectOutliersPCA(miRNASeq1,
#'                             n_components = 2,
#'                             row_threshold = 0.99,
#'                             z_threshold = 3,
#'                             whether_scale = TRUE,
#'                             report_summary = TRUE)
#' result
#' which(result$sample_outliers)
#'
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'   result <- detectOutliersPCA(sample,
#'                               n_components = 2,
#'                               row_threshold = 0.99,
#'                               z_threshold = 3,
#'                               whether_scale = TRUE,
#'                               report_summary = TRUE)
#'   which(result$sample_outliers)
#'   result$miRNA_outliers[1:10, 1:5]
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
#' CRAN Project (2023). \emph{factoextra: Extract and visualize the results of
#' multivariate data analyses.}
#' \href{https://CRAN.R-project.org/package=factoextra}{Link}.
#'
#' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Jolliffe, I. T. (2002). Principal Component Analysis. Springer.
#'
#' @export
#' @importFrom stats prcomp cov mahalanobis
#'
detectOutliersPCA <- function(miRNAdata,
                              n_components = 2,
                              row_threshold = 0.99,
                              z_threshold = 3,
                              whether_scale = TRUE,
                              report_summary = TRUE){

  # Input check
  if (is.matrix(miRNAdata)) {
    miRNAdata <- as.data.frame(miRNAdata)
  } else if (is.data.frame(miRNAdata)) {
    class(miRNAdata) <- "data.frame"
  } else {
    stop("`miRNAdata` must be a data frame or matrix.")
  }

  # Converse non-numeric data
  miRNAdata <- as.data.frame(sapply(miRNAdata, as.numeric))

  # Row-level outliers: PCA + Mahalanobis
  # tag constant columns to avoid error
  constant_cols <- sapply(miRNAdata, function(x) sd(x, na.rm=TRUE) == 0)

  pca_res <- stats::prcomp(miRNAdata[, !constant_cols], scale. = whether_scale)
  scores <- pca_res$x[, 1:n_components, drop = FALSE]

  cov_matrix <- stats::cov(scores)
  center <- colMeans(scores)
  distances <- stats::mahalanobis(scores, center, cov_matrix)
  cutoff <- stats::quantile(distances, row_threshold)
  sample_outliers <- distances > cutoff

  # Column-level outliers: Z-score
  z_scores <- scale(miRNAdata)
  miRNA_outliers <- abs(z_scores) > z_threshold

  # Optional summary
  if (report_summary) {
    message(sprintf("Detected %d row-level outliers (%.1f%% of samples)",
                    sum(sample_outliers), 100*mean(sample_outliers)))
    message(sprintf("Detected %d individual miRNA outliers (%.1f%% of values)",
                    sum(miRNA_outliers), 100*mean(miRNA_outliers)))
  }

  return(list(
    sample_outliers = sample_outliers,
    miRNA_outliers = miRNA_outliers,
    pca_res = pca_res,
    scores = scores
  ))

}


#' Handle Missing Values in miRNA Dataset
#'
#' Detects and imputes missing values in a miRNA dataset. Supports adaptive
#' methods like mean, median, or KNN imputation (Campesato, 2023), ensuring
#' completeness for downstream analysis.
#'
#' @param miRNAdata A numeric data frame or matrix, samples in rows, miRNAs in columns.
#' @param method Method for imputation: "mean", "median", or "knn". Default "median".
#' @param k Number of neighbors for KNN imputation (if method = "knn"). Default 5.
#' @param report_summary Logical, whether to print a summary of missing values before and after imputation. Default TRUE.
#'
#' @return A data frame or matrix with missing values imputed.
#'
#' @references
#' BioConductor Project (2025). \emph{impute: Imputation for microarray data.}
#' \href{https://bioconductor.org/packages/impute/}{Link}.
#'
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Campesato, L. F. (2023). Adaptive missing value imputation methods for omics
#' data. \emph{Bioinformatics}, 39(7), btab789.
#'
#' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Hastie, T., Tibshirani, R., & Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' @export
#' @import stats impute
#'
missingValueHandling <- function(miRNAdata,
                                 method = c("median", "mean", "knn"),
                                 k = 5,
                                 report_summary = TRUE){

  # Check if method selected is valid
  method <- match.arg(method)

  # Input check
  if (is.matrix(miRNAdata)) {
    miRNAdata <- as.data.frame(miRNAdata)
  } else if (is.data.frame(miRNAdata)) {
    class(miRNAdata) <- "data.frame"
  } else {
    stop("`miRNAdata` must be a data frame or matrix.")
  }

  # Converse non-numeric data
  miRNAdata <- as.data.frame(sapply(miRNAdata, as.numeric))

  # Report missing values before
  if (report_summary) {
    message(sprintf("Missing values before imputation: %d",
                    sum(is.na(miRNAdata))))
  }

  # Imputation
  if (method %in% c("mean", "median")) {
    FUN <- ifelse(method == "mean", mean, median)
    miRNAdata <- data.frame(lapply(miRNAdata, function(col) {
      na_idx <- is.na(col)
      if (any(na_idx)) col[na_idx] <- FUN(col, na.rm = TRUE)
      return(col)
    }))
  } else if (method == "knn") {
    # KNN imputation using impute::impute.knn
    if (!requireNamespace("impute", quietly = TRUE)) {
      stop("Package 'impute' is required for KNN imputation.
           Please install it from Bioconductor.")
    }
    # Convert to matrix
    miRNAdata <- as.matrix(miRNAdata)
    imputed <- impute::impute.knn(miRNAdata, k = k)
    miRNAdata <- imputed$data
  }

  # Report missing values after
  if (report_summary) {
    message(sprintf("Missing values after imputation: %d", sum(is.na(miRNAdata))))
  }

  return(miRNAdata)
}
