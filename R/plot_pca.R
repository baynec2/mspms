#' plot_pca
#'
#' Easily create a PCA plot from a QFeatures object containing mspms data. 
#' Ellipses are drawn around the points at a 95% confidence interval.
#' Shape and colors are user specified.
#'
#' @param processed_qf  a QFeatures object containing a SummarizedExperiment
#' named "peptides_norm"
#' @param color the name of the variable you would like to color by.
#' @param shape the name of the variable that you would like to determine 
#' shape by.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' plot_pca(mspms::processed_qf)
plot_pca <- function(processed_qf, color = "time", shape = "condition") {
  mspms_data = mspms_tidy(processed_qf,"peptides_norm")
  # dealing with no visible binding for global variable ‘.’ NOTE
  . <- NULL
  PCA_data <- mspms_data %>%
    dplyr::select("quantCols", "peptide", "group", "condition", "time", "peptides_norm") %>%
    tidyr::pivot_wider(names_from = "peptide", values_from = "peptides_norm", values_fn = NULL) %>%
    # if a peptide has an na remove it
    dplyr::select_if(~ !any(is.na(.)))


  dat <- PCA_data %>%
    dplyr::select(5:length(.))

  # prcomp
  prcomp <- prcomp(dat, center = TRUE, scale = TRUE)


  # Getting PCA values
  PCA_df <- predict(prcomp, dat) %>%
    tibble::as_tibble()

  # Extracting the metadata.
  md <- PCA_data %>%
    dplyr::select(seq_len(4))


  Prop_of_var <- data.frame(summary(prcomp)$importance)[2, ]

  time <- as.factor(md$time)

  condition <- as.factor(md$condition)

  # plotting the PCA
  plot <- PCA_df %>%
    ggplot2::ggplot(ggplot2::aes_string("PC1", "PC2", color = color, linetype = "condition", shape = shape)) +
    ggplot2::geom_point() +
    ggplot2::stat_ellipse(level = 0.95) +
    ggplot2::xlab(paste0("PCA-1 ", round(Prop_of_var[1] * 100, 2), "%")) +
    ggplot2::ylab(paste0("PCA-2 ", round(Prop_of_var[2] * 100, 2), "%"))

  return(plot)
}
