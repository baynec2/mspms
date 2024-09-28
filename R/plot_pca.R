#' plot_pca
#'
#' Easily create a PCA plot from a QFeatures object containing mspms data.
#' Ellipses are drawn around the points at a 95% confidence interval.
#' Shape and colors are user specified.
#'
#' @param mspms_tidy_data tidy mspms data (prepared from QFeatures object
#'  by mspms_tidy)
#' @param value_colname the name of the column containing values.
#' @param color the name of the variable you would like to color by.
#' @param shape the name of the variable that you would like to determine
#' shape by.
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' plot_pca(mspms::mspms_tidy_data)
plot_pca <- function(mspms_tidy_data,
                     value_colname = "peptides_norm",
                     color = "time",
                     shape = "condition") {
    color <- dplyr::sym(color)
    shape <- dplyr::sym(shape)
    all <- prepare_for_PCA(mspms_tidy_data)
    plot <- all[[1]] %>%
        ggplot2::ggplot(ggplot2::aes(.data$PC1, .data$PC2,
            color = !!color,
            linetype = !!shape, shape = !!shape
        )) +
        ggplot2::geom_point() +
        ggplot2::stat_ellipse(level = 0.95) +
        ggplot2::xlab(paste0("PCA-1 ", round(all[[2]][1] * 100, 2), "%")) +
        ggplot2::ylab(paste0("PCA-2 ", round(all[[2]][2] * 100, 2), "%"))

    return(plot)
}
