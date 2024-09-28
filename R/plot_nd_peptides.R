#' plot_nd_peptides
#'
#' plot the percentage of samples each peptide from library was undetected in
#' (if the percentage is > 0).
#'
#' @param processed_qf a QFeatures object containing a SummarizedExperiment
#' named "peptides"
#' @param peptide_library_ids a vector of all peptide library ids in the
#' experiment.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' plot_nd_peptides(mspms::processed_qf)
plot_nd_peptides <- function(processed_qf,
                             peptide_library_ids =
                                 mspms::peptide_library$library_id) {
    nd_peptides <- calc_per_samples_library_nd(
        processed_qf,
        peptide_library_ids
    ) %>%
        dplyr::filter(.data$per_samples_undetected > 0)
    # Plotting
    p1 <- nd_peptides %>%
        ggplot2::ggplot(ggplot2::aes(
            x = stats::reorder(
                .data$library_id,
                .data$per_samples_undetected
            ),
            y = .data$per_samples_undetected
        )) +
        ggplot2::geom_point() +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~ .data$peptide_type, scales = "free_y") +
        ggplot2::xlab("Library ID")
    return(p1)
}
