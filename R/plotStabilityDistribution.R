#' Visualize miRNA Stability Distribution
#'
#' Visualizes the distribution of miRNA stability metrics across all features.
#' Highlights highly stable or highly variable miRNAs, helping users quickly
#' assess feature quality.
#'
#' @param stability_results A list, output from `miRNAStability()`, containing at least `stability_scores` (data.frame with columns: miRNA, CV, MAD)
#' @param metric Character, which stability metric to plot. One of "CV" or "MAD".
#' @param num_top Integer, number of most/least stable miRNAs to highlight. Default 5.
#' @param show_labels Logical, whether to label the highlighted miRNAs on the plot. Default TRUE.
#' @param ... Additional parameters passed to ggplot2 functions.
#'
#' @return A ggplot object visualizing the distribution of stability metrics.
#'
#' @examples
#' # Example 1:
#' # Using miRNASeq1 available with package
#' data(miRNASeq1)
#' result1 <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)
#' plot1 <- plotStabilityDistribution(stability_results = result1)
#'
#' # print(plot1)
#'
#' #Example 2:
#' # Using miRNASeq2 available with package
#' data(miRNASeq2)
#' filled2 <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
#'                                 report_summary = FALSE)
#' result2 <- miRNAStability(miRNAdata = filled2, metrics = "CV",
#'                           report_summary = FALSE)
#' plot2 <- plotStabilityDistribution(stability_results = result2)
#'
#' # print(plot2)
#'
#' \dontrun{
#' # Example 3:
#' # Obtain an external sample miRNASeq dataset (Chodor, 2025)
#' # Example requires the RTCGA.miRNASeq package:
#' if (requireNamespace("RTCGA.miRNASeq", quietly = TRUE)) {
#'   library(RTCGA.miRNASeq)
#'   dim(ACC.miRNASeq) # 240 1048
#'
#'   subset <- RTCGA.miRNASeq::ACC.miRNASeq[1:200, 1:20]
#'   sample <- subset[subset$miRNA_ID == "read_count", ]
#'   result <- miRNAStability(miRNAdata = sample, report_summary = FALSE)
#'   plot <- plotStabilityDistribution(stability_results = result)
#'
#'}}
#'
#' @references
#' Bolstad B. (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Chodor W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Hima Bindu S., et al. (2019). Identification of stable microRNA reference
#' genes for normalization in miRNA expression studies. \emph{BMC Molecular
#' Biology}, 20, 7.
#'
#' Li R., et al. (2016). Evaluation of microRNA expression stability for
#' normalization in quantitative PCR assays. \emph{Biochemical and Biophysical
#' Research Communications}, 471(4), 567-574.
#'
#' Wickham, H. (2016). `ggplot2`: Elegant graphics for data analysis.
#' Springer-Verlag New York. \href{https://ggplot2.tidyverse.org}{Link}.
#'
#' @export
#' @import ggplot2 ggrepel
#' @importFrom utils tail
#'
plotStabilityDistribution <- function(stability_results,
                                      metric = "CV",
                                      num_top = 5,
                                      show_labels = TRUE, ...) {

  # Stop when input is NULL
  if (is.null(stability_results)) stop("Input is NULL")

  # Check if input is a matrix or dataframe
  if (!is.list(stability_results)) {
    stop("`stability_results` must be a list.")
  }

  # Stop when stability_results input is an empty list
  if (length(stability_results) == 0) {
    stop("Empty list input")
  }

  # Check if metric input is valid
  if (!metric %in% c("CV", "MAD")) {
    stop("`metric` must be either 'CV' or 'MAD'.")
  }

  if(is.null(stability_results$stability_scores)){
    stop("Wrong format for input list.")
  }

  scores_df <- stability_results$stability_scores
  if (!metric %in% colnames(scores_df)) {
    stop(paste("Metric", metric, "not found in stability_scores"))
  }

  # Filter miRNA with NA in stability scores (Li et al., 2016)
  scores_df <- scores_df[!is.na(scores_df$CV), ]

  # Stop when scores_df is an empty dataframe after filtering
  if (nrow(scores_df) == 0 || ncol(scores_df) == 0) {
    stop("scores_df is empty")
  }

  # Highlight top/bottom features
  scores_df$highlight <- "Normal"
  sorted_idx <- order(scores_df[[metric]], decreasing = FALSE)
  top_stable <- sorted_idx[1:num_top]
  bottom_stable <- utils::tail(sorted_idx, num_top)
  scores_df$highlight[top_stable] <- "Most Stable"
  scores_df$highlight[bottom_stable] <- "Least Stable"

  # Generate plot using ggplot2 (Wickham, 2016)
  p <- ggplot2::ggplot(scores_df, aes(x = .data[[metric]],
                                      y = reorder(scores_df$miRNA, .data[[metric]]),
                                      color = scores_df$highlight)) +
                                      ggplot2::geom_point(size = 3) +
    scale_color_manual(values = c("Most Stable" = "blue",
                                  "Least Stable" = "red",
                                  "Normal" = "gray")) +
    labs(x = paste(metric, "across miRNAs"),
         y = "miRNA",
         color = "Feature stability") +
    theme_minimal()

  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_text_repel(data = subset(scores_df,
                                                    scores_df$highlight != "Normal"),
                                      aes(label = scores_df$miRNA),
                                      size = 3,
                                      max.overlaps = Inf)
  }

  return(p)

}

