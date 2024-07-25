#' plot_nd_peptides
#'
#' plot peptides from libray that were not detected in the data
#'
#' @param nd_peptides a tibble containing info on peptides not detected by the 
#' find_nd_peptides function.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' nd_peptides <- find_nd_peptides(
#'   mspms::peaks_prepared_data,
#'   mspms::peptide_library,
#'   mspms::design_matrix
#' )
#' plot_nd_peptides(nd_peptides)
plot_nd_peptides <- function(nd_peptides) {
  # Calculating number of total samples in study
  n_samples <- length(unique(nd_peptides$sample))

  # Summarizing per peptide
  sum <- nd_peptides %>%
    dplyr::group_by(.data$missing_peptide_id,.data$type_missing) %>%
    dplyr::summarise(
      n_missing = n(),
      percent_samples_undetected =
        round((.data$n_missing / n_samples * 100), 1)
    )
  # Plotting
  p1 <- sum %>%
    ggplot2::ggplot(ggplot2::aes(
      x = stats::reorder(
        .data$missing_peptide_id,
        .data$percent_samples_undetected
      ),
      y = .data$percent_samples_undetected
    )) +
    ggplot2::geom_point() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~.data$type_missing, scales = "free_y") +
    ggplot2::xlab("Library ID")
  return(p1)
}
