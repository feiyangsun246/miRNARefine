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
#'
#' @return A list where each element is the normalized dataset corresponding to a method.
#'
#' @examples
#' data(toy_miRNA)
#' normalized <- compareNormalization(toy_miRNA)
#' normalized$log2[1:5, 1:5]
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
#' @importFrom stats scale
