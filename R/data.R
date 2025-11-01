#' miRNA Sequence dataset from TCGA project
#'
#' A small subset of the BRCA.miRNASeq dataset from the RTCGA.miRNASeq package.
#' This dataset contains expression information for 10 samples and 18 miRNAs.
#' Each value in the data frame represents **reads per million miRNA mapped (RPM)** for that sample and miRNA.
#' It is provided as an example dataset for demonstration and testing purposes.
#'
#' @source \url{https://bioconductor.org/packages/RTCGA.miRNASeq/} (RTCGA.miRNASeq package)
#'
#' @format A data frame with 18 columns and 10 rows:
#' \describe{
#'   \item{miRNA1}{RPM value for miRNA1 in each sample}
#'   \item{miRNA2}{RPM value for miRNA2 in each sample}
#'   ...
#'   \item{miRNA18}{RPM value for miRNA18 in each sample}
#' }
#'
#' @examples
#' \dontrun{
#' miRNASeq1
#' }
"miRNASeq1"

#' Toy miRNA Dataset
#'
#' A small toy miRNA dataset including missing values and outliers.
#' This dataset is designed for testing preprocessing and filtering functions.
#' Each row represents a sample, and each column represents a miRNA feature.
#' Values represent **reads per million miRNA mapped (RPM)**.
#' Some values are intentionally set as \code{NA} or extreme outliers
#' to simulate noisy, unclean data.
#'
#' @format A data frame with 20 columns and 30 rows:
#' \describe{
#'   \item{miRNA1}{RPM value for miRNA1 in each sample}
#'   \item{miRNA2}{RPM value for miRNA2 in each sample}
#'   ...
#'   \item{miRNA20}{RPM value for miRNA20 in each sample}
#' }
#'
#' @examples
#' \dontrun{
#' miRNASeq2
#' }
"miRNASeq2"
