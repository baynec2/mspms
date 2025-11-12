#' plot_volcano
#'
#' create a volcano plot to generate log2fc and adjusted p values for
#' experimental conditions
#'
#' @param log2fc_t_test_data a tibble containing the log2fc
#' and adjusted p values
#' @param log2fc_threshold  the log2fc threshold that you want displayed on plot
#' @param padj_threshold  the padj threshold that you want displayed on plot
#' @param facets how facets should be displayed. Accepted values are grid and
#' wrap
#' @param ncol ncol to include if facets = "wrap"
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' p1 <- mspms::plot_volcano(mspms::log2fc_t_test_data, log2fc_threshold = 3)
#' p1
plot_volcano <- function(log2fc_t_test_data,
                         log2fc_threshold = 3,
                         padj_threshold = 0.05,
                         facets = "grid",
                         ncol = 1) {
  ncol <- length(unique(log2fc_t_test_data$time))
  p1 <- log2fc_t_test_data %>%
    dplyr::mutate(`log_10_p.adj` = -log10(.data$p.adj)) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$log2fc,
      y = .data$log_10_p.adj
    )) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = "dashed", color = "red"
    ) +
    ggplot2::geom_vline(
      xintercept = log2fc_threshold,
      linetype = "dashed", color = "red"
    ) +
    ggplot2::geom_vline(
      xintercept = -log2fc_threshold,
      linetype = "dashed", color = "red"
    ) +
    ggplot2::labs(x = "Log2 (Time TX/T0)", y = "-log10(p.adj value)")
  if (facets == "grid") {
    p2 <- p1 + ggplot2::facet_grid(
      rows = dplyr::vars(.data$condition),
      cols = dplyr::vars(.data$time)
    )
  } else if (facets == "wrap") {
    p2 <- p1 + ggplot2::facet_wrap(~ .data$condition + .data$time,
      ncol = ncol
    )
  } else {
    stop("facets must be equal to either grid or wrap")
  }
  return(p2)
}
