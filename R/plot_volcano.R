#' plot_volcano
#'
#' create a volcano plot to generate log2fc and adjusted p values for
#' experimental conditions
#'
#' @param log2fc_t_test_data = this is the data frame that contains the log2fc
#' and adjusted p values
#' @param log2fc_threshold  = this is the threshold that you want displayed
#'  on plot
#' @param facets = grid or wrap. If grid facets for each condition and time are
#' shown even if data does not exist. Wrap on the other hand is condensed. 
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Generating log2fc t test data
#' log2fc_t_test_data <- mspms::log2fc_t_test(mspms::mspms_data)
#' # plotting data.
#' p1 <- mspms::plot_volcano(log2fc_t_test_data, log2fc_threshold = 3)
#' p1
plot_volcano <- function(log2fc_t_test_data,
                         log2fc_threshold = 3,
                         facets = "grid") {
  ncol <- length(unique(log2fc_t_test_data$time))
  p1 <- log2fc_t_test_data %>%
    dplyr::mutate(`log_10_p.adj` = -log10(.data$p.adj)) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "log2fc", y = "log_10_p.adj")) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
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
  if(facets == "grid"){
    p2 = p1 + ggplot2::facet_grid(rows = vars(condition),
                                  cols = vars(time)) 
  } else if (facets == "wrap"){
    p2 = p1 + ggplot2::facet_wrap(~condition + time,
                                  ncol = ncol) 
  } else {
    stop("facets must be either grid or wrap")
  }
  return(p2)
}
