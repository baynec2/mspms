#' plot_cleavages_per_pos
#'
#' make a plot of the count of cleavages per position.
#'
#' @param count_of_cleavages = tibble with number of cleavages per position, for
#' each condition and time point
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'c = data.frame(time = rep(60,13),
#'                   condition = rep("DMSO",13),
#'                   cleavage_pos = 1:13,
#'                   n = 1:13)
#' plot_cleavages_per_pos(c)
plot_cleavages_per_pos = function(count_of_cleavages){

  p1 = ggplot2::ggplot(count_of_cleavages,
                       ggplot2::aes_string(x = "cleavage_pos",
                                           y = "n")) +
    ggplot2::geom_point() +
    ggplot2::geom_line()+
    ggplot2::facet_wrap(~condition + time, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(0, 13))

  return(p1)
}


