#' plot_cleavages_per_pos
#'
#' plot the number of cleavages at each
#'
#' @param sig_cleavage_data a tibble of data of interest containing a column
#' labeled peptide, cleavage_seq, condition, and cleavage_pos.
#' @param ncol the number of columns to plot.
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Defining the significant peptides
#' sig_cleavage_data <- log2fc_t_test_data %>%
#'   dplyr::filter(p.adj <= 0.05, log2fc > 3)
#' # Plotting
#' p1 <- mspms::plot_cleavages_per_pos(sig_cleavage_data)
#' p1
plot_cleavages_per_pos <- function(sig_cleavage_data,
                                   ncol = NULL) {
  count_cleavages_per_pos <- count_cleavages_per_pos(sig_cleavage_data)
  p1 <- count_cleavages_per_pos %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$cleavage_pos,
      y = .data$n,
      color = .data$time
    )) +
    ggplot2::facet_wrap(~ .data$condition, scales = "free_y", ncol = ncol) +
    ggplot2::geom_point() +
    ggplot2::geom_line()
  return(p1)
}
