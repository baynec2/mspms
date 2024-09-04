#' plot_rt_qc
#' plot the retention times of detected peptides as a histogram or density plot
#' @param prepared_data  = data prepared by mspms package
#' @param design_matrix  = design matrix
#' @param type = type of plot - histogram or density
#' @param binwidth = binwidth for histogram
#' @param facet_by_group = whether to facet by group, TRUE or FALSE
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'plot_rt_qc(mspms::peaks_prepared_data,
#'           mspms::design_matrix,
#'           type = "density",
#'           facet_by_group = TRUE)
plot_rt_qc <- function(prepared_data,
                       design_matrix,
                       type = "histogram",
                       binwidth = 2,
                       facet_by_group = TRUE) {
  # Plot the QC of the RTs
  long <- prepared_data %>%
    tidyr::pivot_longer(4:length(.), names_to = "sample") %>%
    dplyr::filter(!is.na(value)) %>% 
    dplyr::inner_join(design_matrix, by = "sample")
  
  # Formatting data in the long format
  p1 <- ggplot2::ggplot(long, ggplot2::aes(x = RT)) +
    if (type == "histogram") {
      ggplot2::geom_histogram(binwidth = binwidth) 
    } else if (type == "density") {
      ggplot2::geom_density()
    } else {
      stop("type must be histogram or density")
    }
  
  if(isTRUE(facet_by_group)) {
    p1 <- p1 + ggplot2::facet_wrap(~.data$group, scales = "free_y")
  }else if(isFALSE(facet_by_group)) {
    p1 <- p1
  }else {
    stop("facet_by_group must be TRUE or FALSE")
  }

  return(p1)
}