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
#'   \item{hsa-let-7a-1}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7a-2}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7a-3}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7b}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7c}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7d}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7e}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7f-1}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7f-2}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7g}{RPM value for the miRNA in each sample}
#'   \item{hsa-let-7i}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-1-1}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-1-2}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-100}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-101-1}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-101-2}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-103-1}{RPM value for the miRNA in each sample}
#'   \item{hsa-mir-103-1-as}{RPM value for the miRNA in each sample}
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
#'   \item{miRNA3}{RPM value for miRNA3 in each sample}
#'   \item{miRNA4}{RPM value for miRNA4 in each sample}
#'   \item{miRNA5}{RPM value for miRNA5 in each sample}
#'   \item{miRNA6}{RPM value for miRNA6 in each sample}
#'   \item{miRNA7}{RPM value for miRNA7 in each sample}
#'   \item{miRNA8}{RPM value for miRNA8 in each sample}
#'   \item{miRNA9}{RPM value for miRNA9 in each sample}
#'   \item{miRNA10}{RPM value for miRNA10 in each sample}
#'   \item{miRNA11}{RPM value for miRNA11 in each sample}
#'   \item{miRNA12}{RPM value for miRNA12 in each sample}
#'   \item{miRNA13}{RPM value for miRNA13 in each sample}
#'   \item{miRNA14}{RPM value for miRNA14 in each sample}
#'   \item{miRNA15}{RPM value for miRNA15 in each sample}
#'   \item{miRNA16}{RPM value for miRNA16 in each sample}
#'   \item{miRNA17}{RPM value for miRNA17 in each sample}
#'   \item{miRNA18}{RPM value for miRNA18 in each sample}
#'   \item{miRNA19}{RPM value for miRNA19 in each sample}
#'   \item{miRNA20}{RPM value for miRNA20 in each sample}
#' }
#'
#' @examples
#' \dontrun{
#' miRNASeq2
#' }
"miRNASeq2"
