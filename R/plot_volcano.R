#' plot_volcano
#'
#' create a volcano plot to generate log2fc and adjusted p values for experimental conditions
#'
#' @param log2fc_t_test_data = this is the data frame that contains the log2fc and adjusted p values
#' @param log2fc_threshold  = this is the threshold that you want displayed on plot
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Generating log2fc t test data
#' log2fc_t_test_data = mspms::log2fc_t_test(mspms::prepared_for_stats)
#' # plotting data.
#'p1 =  mspms::plot_volcano(log2fc_t_test_data,log2fc_threshold = 3)
#'p1
plot_volcano = function(log2fc_t_test_data,log2fc_threshold = 3){

  ncol = length(unique(log2fc_t_test_data$group2))

  p1 = log2fc_t_test_data %>%
    ggplot2::ggplot(ggplot2::aes(x = log2fc,y = -log10(p.adj)))+
    ggplot2::geom_point(size =0.5)+
    ggplot2::geom_hline(yintercept = -log10(0.05),linetype = "dashed",color = "red")+
    ggplot2::geom_vline(xintercept = log2fc_threshold, linetype = "dashed",color = "red")+
    ggplot2::geom_vline(xintercept = -log2fc_threshold, linetype = "dashed",color = "red")+
    ggplot2::labs(x = "Log2 (Time TX/T0)",y = "-log10(p value)")+
    ggplot2::facet_wrap(~ condition + group2,ncol = ncol)

  return(p1)
}
