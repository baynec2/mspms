#' plot_qc_check
#' plot the the percentage of the peptide library undetected in each sample per
#' each sample group.
#' @param processed_qf QFeatures object containing a SummarizedExperiment
#' named "peptides"
#' @param peptide_library a vector of all peptide library ids in the
#' experiment.
#' @param full_length_threshold percent to use as threshold visualized as a
#' vertical blue dashed line
#' @param cleavage_product_threshold percent to use as a threshold visualized as
#' a red dashed line
#' @param ncol n columns.
#' @return a ggplot2 object.
#' @export
#'
#' @examples
#' plot_qc_check(mspms::processed_qf)
plot_qc_check <- function(processed_qf,
                          peptide_library = mspms::peptide_library$library_id,
                          full_length_threshold = NULL,
                          cleavage_product_threshold = NULL,
                          ncol = 2) {
  if (!is.character(peptide_library)) {
    stop("peptide library musyt be a charachter vector")
  }
  qc_check_data <- prepare_qc_check_data(
    processed_qf,
    peptide_library
  )
  # Plotting per
  p1 <- qc_check_data %>%
    ggplot2::ggplot(ggplot2::aes(
      color = .data$peptide_type,
      fill = .data$peptide_type
    )) +
    ggplot2::geom_histogram(ggplot2::aes(.data$per_library_id_undetected),
      binwidth = 0.1,
      alpha = 0.5
    ) +
    ggplot2::geom_density(ggplot2::aes(.data$per_library_id_undetected),
      alpha = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = full_length_threshold, linetype = "dashed",
      color = "#00BFC4"
    ) +
    ggplot2::geom_vline(
      xintercept = cleavage_product_threshold, linetype = "dashed",
      color = "#F8766D"
    ) +
    ggplot2::facet_wrap(~ .data$group, scales = "free_y", ncol = ncol) +
    ggplot2::ggtitle("Percentage of Library ID Undetected In Each Sample") +
    ggplot2::ylab("count")
  return(p1)
}
