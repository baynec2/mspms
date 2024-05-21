#' plot_pca
#'
#' this creates a pca plot of the data where the points are colored by time and shaped by condition. elipses are drawn around the points at a 95% confidence interval.
#'
#' @param prepared_for_stats = this is the data that has been prepared for stats by the prepare_for_stats() function
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#'
#'plot_pca(mspms::prepared_for_stats)
#'
#'
plot_pca = function(prepared_for_stats){
  # dealing with no visible binding for global variable ‘.’ NOTE
  . = NULL
  PCA_data = prepared_for_stats %>%
    dplyr::select(.data$sample,.data$Peptide,.data$group,.data$condition,.data$time,.data$value) %>%
    tidyr::pivot_wider(names_from = .data$Peptide,values_from = .data$value,values_fn = mean) %>%
    # if a peptide has an na remove it
    dplyr::select_if(~ !any(is.na(.)))


  dat = PCA_data %>%
    dplyr::select(5:length(.))

  # prcomp
  prcomp = prcomp(dat, center = TRUE, scale = TRUE)


  #Getting PCA values
  PCA_df = predict(prcomp,dat) %>%
    tibble::as_tibble()

  # Extracting the metadata.
  md = PCA_data %>%
    dplyr::select(1:4)


  Prop_of_var = data.frame(summary(prcomp)$importance)[2,]

  time = as.factor(md$time)

  condition = md$condition

  #plotting the PCA
  plot = PCA_df %>%
    ggplot2::ggplot(ggplot2::aes_string("PC1","PC2",color = "time", linetype = "condition",shape = "condition"))+
    ggplot2::geom_point()+
    ggplot2::stat_ellipse(level=0.95)+
    ggplot2::xlab(paste0("PCA-1 ",round(Prop_of_var[1]*100,2),"%"))+
    ggplot2::ylab(paste0("PCA-2 ",round(Prop_of_var[2]*100,2),"%"))+
    ggplot2::theme_bw()

  return(plot)

}
