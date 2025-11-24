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
#' # Obtain an external sample miRNASeq dataset (Chodor, 2025)
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   subset <- RTCGA.miRNASeq::ACC.miRNASeq[1:200, 1:20]
#'   sample_ori <- subset[subset$miRNA_ID == "read_count", ]
#'
#'   # deal with missing values first in case of errors
#'   sample <- missingValueHandling(miRNAdata = sample_ori,
#'                                  method = "knn",
#'                                  k = 5,
#'                                  report_summary = TRUE))
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
#' R Core Team. (2025). Package `stats`. R: A Language and Environment for
#' Statistical Computing. R Foundation for Statistical Computing, Vienna,
#' Austria <https://www.R-project.org/>
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
  miRNAdata <- inputCheckGeneral(miRNAdata)

  # Stop when there is NA in the dataset
  if (any(is.na(miRNAdata))) {
    stop("Dataset contains missing values (NA).
         Consider running missingValueHandling() first.")
  }

  # Stop when n_components is not valid
  if (!is.numeric(n_components) || length(n_components) != 1 || n_components <= 0) {
    stop("`n_components` must be a single positive number.")
  }

  # Stop when row_threshold is not valid
  if (!is.numeric(row_threshold) || length(row_threshold) != 1 ||
      row_threshold <= 0 || row_threshold > 1) {
    stop("`row_threshold` must be a single number in (0, 1].")
  }

  # Stop when z_threshold is not valid
  if (!is.numeric(z_threshold) || length(z_threshold) != 1 || z_threshold <= 0) {
    stop("`z_threshold` must be a single positive number.")
  }

  # Row-level outliers: PCA + Mahalanobis (Jolliffe, 2002; CRAN Project, 2023)
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

  # Perform outliers calculation (R Core Team, 2025)
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
#' # Obtain an external sample miRNASeq dataset (Chodor, 2025)
#' # Example requires the RTCGA.miRNASeq package:
#' \dontrun{
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   subset <- RTCGA.miRNASeq::ACC.miRNASeq[1:200, 1:20]
#'   sample <- subset[subset$miRNA_ID == "read_count", ]
#'   filled <- missingValueHandling(miRNAdata = sample,
#'                                  method = "knn",
#'                                  k = 5,
#'                                  report_summary = TRUE)
#'   head(filtered)
#' }}
#'
#' @references
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
#' Hastie, T., Tibshirani, R., and Friedman, J. (2009).
#' \emph{The Elements of Statistical Learning}, 2nd Edition. Springer.
#'
#' Hastie T., Tibshirani R., Narasimhan B., and Chu G. (2025). impute: impute:
#' Imputation for microarray data. R package version 1.84.0,
#' \href{https://bioconductor.org/packages/impute/}{Link}.
#'
#' Robinson, M. D., McCarthy, D. J., and Smyth, G. K. (2010). edgeR: a
#' Bioconductor package for differential expression analysis of digital gene
#' expression data. \emph{Bioinformatics}, 26, 139-140.
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/}{Link}.
#'
#' @export
#' @import stats impute
#'
missingValueHandling <- function(miRNAdata,
                                 method = "median",
                                 k = 5,
                                 report_summary = TRUE){

  # Check if method is valid
  if (!method %in% c("mean", "median", "knn")) {
    stop("`method` must be one of 'mean', 'median', or 'knn'.")
  }

  # Check if input is a matrix or dataframe
  miRNAdata <- inputCheckGeneral(miRNAdata)

  # Check if any column is entirely NA
  if (any(colSums(is.na(miRNAdata)) == nrow(miRNAdata))) {
    stop("One or more columns contain only NA values.
         Consider remove such columns first.")
  }

  # Report missing values before
  if (report_summary) {
    message(sprintf("Missing values before imputation: %d",
                    sum(is.na(miRNAdata))))
  }

  # Imputation using mean or median methods (Hastie et al., 2009; Campesato, 2023)
  if (method %in% c("mean", "median")) {
    FUN <- ifelse(method == "mean", mean, median)
    filled_miRNA <- data.frame(lapply(miRNAdata, function(col) {
      na_idx <- is.na(col)
      if (any(na_idx)) col[na_idx] <- FUN(col, na.rm = TRUE)
      return(col)
    }))

  } else if (method == "knn") {
    # KNN imputation using impute::impute.knn (Hastie et al., 2025)
    if (!requireNamespace("impute", quietly = TRUE)) {
      stop("Package 'impute' is required for KNN imputation.
           Please install it from Bioconductor.")
    }

    if (!is.numeric(k) || k <= 0 || k != round(k)) {
      stop("`k` must be a positive integer for kNN imputation.")
    }

    if (k >= nrow(miRNAdata)) {
      stop("`k` must be smaller than the number of samples.")
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

