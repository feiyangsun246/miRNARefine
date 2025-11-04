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
#' @details
#' Before running `detectOutliersPCA()`, ensure that the dataset contains
#' no missing values.
#' You may use `missingValueHandling()` to fill missing data to avoid errors.
#'
#' @examples
#' # Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq1)
#' result <- detectOutliersPCA(miRNAdata = miRNASeq1,
#'                             n_components = 2,
#'                             row_threshold = 0.99,
#'                             z_threshold = 3,
#'                             whether_scale = TRUE,
#'                             report_summary = TRUE)
#' result
#' which(result$sample_outliers)
#'
#' \dontrun{
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample_ori <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'
#'   # deal with missing values first in case of errors
#'   sample <- missingValue
#'
#'   result <- detectOutliersPCA(miRNAdata = sample,
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

  # Row-level outliers: PCA + Mahalanobis
  # tag constant columns to avoid error
  constant_cols <- sapply(miRNAdata, function(x) sd(x, na.rm=TRUE) == 0)

  pca_res <- stats::prcomp(miRNAdata[, !constant_cols], scale. = whether_scale)
  scores <- pca_res$x[, 1:n_components, drop = FALSE]

  cov_matrix <- stats::cov(scores)

  # Check if the covariance matrix is singular or nearly singular
  if (any(is.na(cov_matrix)) || det(cov_matrix) < .Machine$double.eps) {
    warning("Covariance matrix is singular or ill-conditioned;
            skipping Mahalanobis computation.")
    return(data.frame(Sample = rownames(miRNAdata), IsOutlier = FALSE))
  }

  center <- colMeans(scores)
  distances <- stats::mahalanobis(scores, center, cov_matrix)
  cutoff <- stats::quantile(distances, row_threshold)
  sample_outliers <- distances > cutoff

  # Column-level outliers: Z-score
  z_scores <- scale(miRNAdata)
  miRNA_outliers <- abs(z_scores) > z_threshold

  # Optional summary
  if (isTRUE(report_summary)) {
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
#' @examples
#' # Example 1:
#' # Using miRNASeq data available with package
#' data(miRNASeq2)
#' filled <- missingValueHandling(miRNAdata = miRNASeq2,
#'                                method = "median",
#'                                k = 5,
#'                                report_summary = TRUE)
#' head(filled)
#'
#' # Example 2:
#' # Obtain an external sample miRNASeq dataset
#' # Example requires the RTCGA.miRNASeq package:
#' \dontrun{
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   sample <- RTCGA.miRNASeq::ACC.miRNASeq[1:100, 1:20]
#'   filled <- missingValueHandling(miRNAdata = sample,
#'                                  method = "knn",
#'                                  k = 5,
#'                                  report_summary = TRUE)
#'   head(filtered)
#' }}
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
#' Robinson, M. D., McCarthy, D. J. and Smyth, G. K. (2010). edgeR: a
#' Bioconductor package for differential expression analysis of digital gene
#' expression data. \emph{Bioinformatics}, 26, 139-140.
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/}{Link}
#'
#' @export
#' @import stats impute
#'
missingValueHandling <- function(miRNAdata,
                                 method = c("median", "mean", "knn"),
                                 k = 5,
                                 report_summary = TRUE){

  # Match method
  method <- match.arg(method)

  # Check if input is valid
  if (is.matrix(miRNAdata)) {
    miRNAdata <- as.data.frame(miRNAdata)
  } else if (is.data.frame(miRNAdata)) {
    class(miRNAdata) <- "data.frame"
  } else {
    stop("`miRNAdata` must be a data frame or matrix")
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

  # Report missing values before
  if (report_summary) {
    message(sprintf("Missing values before imputation: %d",
                    sum(is.na(miRNAdata))))
  }

  # Imputation
  if (method %in% c("mean", "median")) {
    FUN <- ifelse(method == "mean", mean, median)
    filled_miRNA <- data.frame(lapply(miRNAdata, function(col) {
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
    filled_miRNA <- imputed$data
  }

  # Report missing values after
  if (isTRUE(report_summary)) {
    message(sprintf("Missing values after imputation: %d",
                    sum(is.na(miRNAdata))))
  }

  return(filled_miRNA)

}
