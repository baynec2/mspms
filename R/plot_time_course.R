#' plot_time_course
#'
#' Easily plot a time course of all peptides in a QFeatures object by
#' peptide.
#' 
#'
#' @param processed_qf a QFeatures object containing a SummarizedExperiment
#' named "peptides_norm".
#' @param summarize_by_mean whether to summarise by mean (TRUE- show error bars 
#' +- 1 standard deviation) or not (FALSE)
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Determining peptide of interest
#'max_log2fc_pep = mspms::log2fc_t_test_data %>% 
#' dplyr::filter(p.adj <= 0.05,log2fc > 3) %>% 
#' dplyr::filter(log2fc ==max(log2fc)) %>% 
#' dplyr::pull(peptide)
#' 
#' # Defining QFeatures filter
#' filter = QFeatures::VariableFilter(field = "peptide",value = max_log2fc_pep,
#'condition = "==")
#'
#'# ploting time course
#'plot_time_course(sel,summarize_by_mean = FALSE)
plot_time_course <- function(processed_qf,
                             summarize_by_mean = FALSE) {
  mspms_data = mspms_tidy(processed_qf) 
  
  if(isTRUE(summarize_by_mean)){
  p1 <- mspms_data %>%
    dplyr::group_by(.data$peptide, .data$condition, .data$time,
                    .data$cleavage_seq) %>%
    dplyr::summarize(mean = mean(.data$peptides_norm, na.rm = TRUE),
                     sd = sd(.data$peptides_norm, na.rm = TRUE)) %>%
    ggplot2::ggplot(ggplot2::aes(x =.data$time,
                                 y = .data$mean,
                                 color = .data$condition)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymax = mean + sd,
      ymin = mean - sd
    ), width = 15) +
    ggplot2::facet_wrap(~peptide+cleavage_seq, scales = "free_y") }
  else{
    p1 <- mspms_data %>%
      ggplot2::ggplot(ggplot2::aes(x =.data$time,
                                   y = .data$peptides_norm,
                                   color = .data$condition
                                   )) +
      ggplot2::geom_point() +
      ggplot2::stat_summary(ggplot2::aes(group=.data$condition), fun =mean,
                            geom="line")+
      ggplot2::facet_wrap(~peptide+cleavage_seq, scales = "free_y")
  }

  return(p1)
}
