#' plot_time_course
#'
#' this function plots the time course of each peptide in the data set.
#'
#' @param prepared_for_stats this is a data frame that has been prepared for stats.
#' It should have the following columns: time, value, condition, and Peptide. It is in long format.
#' Can be filtered to only include the peptides you want to plot.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#' p1 = mspms::prepared_for_stats %>%
#' dplyr::filter(Peptide == "A_GLFNYNQLRGF") %>%
#' mspms::plot_time_course()
#'
#' p1

plot_time_course = function(prepared_for_stats){

 p1 = prepared_for_stats %>%
  dplyr::group_by(Peptide,condition,time) %>%
  dplyr::summarize(mean = mean(value,na.rm=T),sd = sd(value,na.rm = T)) %>%
  ggplot2::ggplot(ggplot2::aes(x = time, y = mean, color = condition)) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::geom_errorbar(ggplot2::aes(ymax = mean + sd, ymin = mean - sd),width = 15) +
  ggplot2::facet_wrap(~Peptide, scales = "free_y") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "bottom")


  return(p1)

}



