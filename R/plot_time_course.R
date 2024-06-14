#' plot_time_course
#'
#' this function plots the time course of each peptide in the data set.
#'
#' @param mspms_data this is data from the mspms pipeline.
#' It should have the following columns: time, value, condition, and Peptide. It is in long format.
#' Can be filtered to only include the peptides you want to plot.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' p1 <- mspms::mspms_data %>%
#'   dplyr::filter(Peptide == "A_GLFNYNQLRGF") %>%
#'   mspms::plot_time_course()
#'
#' p1
plot_time_course <- function(mspms_data) {
  p1 <- mspms_data %>%
    dplyr::group_by(.data$Peptide, .data$condition, .data$time,.data$cleavage_seq) %>%
    dplyr::summarize(mean = mean(.data$value, na.rm = TRUE), sd = sd(.data$value, na.rm = TRUE)) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "time", y = "mean", color = "condition")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymax = mean + sd,
      ymin = mean - sd
    ), width = 15) +
    ggplot2::facet_wrap(~Peptide+cleavage_seq, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  return(p1)
}
