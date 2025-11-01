#' miRNA Sequence dataset from TCGA project
#'
#' A small subset of the BRCA.miRNASeq dataset from the RTCGA.miRNASeq package.
#' This dataset contains expression information for 10 samples and 18 miRNAs.
#' Each value in the data frame represents **reads per million miRNA mapped (RPM)** for that sample and miRNA.
#' It is provided as an example dataset for demonstration and testing purposes.
#'
#' @source \url{https://bioconductor.org/packages/RTCGA.miRNASeq/} (RTCGA.miRNASeq package)
#'
#' @format A data frame with 10 rows and 18 columns:
#' \describe{
#'   \item{Sample1}{Expression value for miRNA 1 in sample 1}
#'   \item{Sample2}{Expression value for miRNA 2 in sample 1}
#'   ...
#'   \item{Sample10}{Expression value for miRNA 18 in sample 1}
#' }
#'
#' @examples
#' \dontrun{
#' miRNASeq1
#' }
"miRNASeq1"
