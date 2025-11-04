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
#' @references
#' Bolstad B (2025). preprocessCore: A collection of pre-processing functions.
#' R package version 1.72.0.
#' \href{https://bioconductor.org/packages/preprocessCore}{Link}.
#'
#' Chodor, W. (2025). \emph{RTCGA.miRNASeq: miRNASeq datasets from The Cancer
#' Genome Atlas Project}. R package version 1.36.0.
#'
#' Hima Bindu, S., et al. (2019). Identification of stable microRNA reference
#' genes for normalization in miRNA expression studies. \emph{BMC Molecular
#' Biology}, 20, 7.
#'
#' Li, R., et al. (2016). Evaluation of microRNA expression stability for
#' normalization in quantitative PCR assays. \emph{Biochemical and Biophysical
#' Research Communications}, 471(4), 567-574.
#'
#' @export
#' @import ggplot2
#' @importFrom utils tail
#'
plotStabilityDistribution <- function(stability_results,
                                      metric = "CV",
                                      num_top = 5,
                                      show_labels = TRUE, ...) {

  # Check if metric input is valid
  if (!metric %in% c("CV", "MAD")) {
    stop("`metric` must be either 'CV' or 'MAD'")
  }

  scores_df <- stability_results$stability_scores
  if (!metric %in% colnames(scores_df)) {
    stop(paste("Metric", metric, "not found in stability_scores"))
  }

  # Highlight top/bottom features
  scores_df$highlight <- "Normal"
  sorted_idx <- order(scores_df[[metric]], decreasing = FALSE)
  top_stable <- sorted_idx[1:num_top]
  bottom_stable <- utils::tail(sorted_idx, num_top)
  scores_df$highlight[top_stable] <- "Most Stable"
  scores_df$highlight[bottom_stable] <- "Least Stable"

  p <- ggplot2::ggplot(scores_df, aes(x = .data[[metric]],
                                      y = reorder(miRNA, .data[[metric]]),
                                      color = highlight)) +
                                      geom_point(size = 3) +
    scale_color_manual(values = c("Most Stable" = "blue",
                                  "Least Stable" = "red",
                                  "Normal" = "gray")) +
    labs(x = paste(metric, "across miRNAs"),
         y = "miRNA",
         color = "Feature stability") +
    theme_minimal()

  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_text_repel(data = subset(scores_df, highlight != "Normal"),
                                      aes(label = miRNA),
                                      size = 3,
                                      max.overlaps = Inf)
  }

}

