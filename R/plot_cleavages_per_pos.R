#' plot_cleavages_per_pos
#'
#' make a plot of the count of cleavages per position.
#'
#' @param count_of_cleavages = tibble with number of cleavages per position, for
#' each condition and time point
#' @param facets = condition, grid, or wrap If grid facets for each condition and 
#' time are shown even if data does not exist. 
#' Wrap on the other hand is condensed. None shows plots overlayed by time with
#' different colors.
#' 
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
plot_cleavages_per_pos = function(count_of_cleavages,
                                  facets = "condition"){
  p1 = ggplot2::ggplot(count_of_cleavages,
                       ggplot2::aes_string(x = "cleavage_pos",
                                           y = "n")) +
    ggplot2::geom_point() +
    ggplot2::geom_line()+
    ggplot2::scale_x_continuous(breaks = seq(0, 13))

  if(facets == "grid"){
    p2 = p1 + ggplot2::facet_grid(rows = dplyr::vars(condition),
                                  cols = dplyr::vars(time)) 
  } else if (facets == "wrap"){
    p2 = p1 + ggplot2::facet_wrap(~condition + time,
                                  ncol = ncol,
                                  scales = "free_y") 
  } else if(facets == "condition"){
    p2 = ggplot2::ggplot(count_of_cleavages,
                         ggplot2::aes_string(x = "cleavage_pos",
                                             y = "n",
                                             color = "time"))+
      ggplot2::geom_point()+
      ggplot2::geom_line()+
      ggplot2::scale_x_continuous(breaks = seq(0, 13))+
      ggplot2::facet_wrap(~condition)
  } else {
    stop("facets must be either condition,grid or wrap")
  }
  return(p2)
}



