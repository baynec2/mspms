#' plot_heatmap
#'
#' This produces a heatmaply interactive heatmap of the QFeatures object with
#' color bars representing the condition and time for each sample in each row.
#'
#' Each column has a colored bar representing whether the peptide is a
#' cleavage product or a full length member of the peptide library.
#'
#' @param mspms_tidy_data tidy mspms data (prepared from QFeatures object
#' by mspms_tidy())
#' @param value_colname the name of the column containing values.
#' @param scale how would you like the data scaled? default is none,
#'  but can also be "row", "column", or "none"
#' @param plot_method what plot method would you like to use, can use
#' plotly or ggplot2.
#' @param show_dendrogram Logical vector of length two, controlling whether
#' the row and/or column dendrograms are displayed. If a logical scalar is
#' provided, it is repeated to become a logical vector of length two.

#' @return a heatmaply interactive heatmap
#' @export
#'
#' @examples
#' plot_heatmap(mspms::mspms_tidy_data)
plot_heatmap <- function(mspms_tidy_data,
                         value_colname = "peptides_norm",
                         scale = "column",
                         plot_method = "plotly",
                         show_dendrogram = c(TRUE, TRUE)) {
  value_colname <- dplyr::sym(value_colname)
  heatmap_data <- mspms_tidy_data %>%
    dplyr::select(
      "quantCols", "condition", "time",
      "peptide", !!value_colname
    ) %>%
    tidyr::pivot_wider(
      names_from = "peptide",
      values_from = "peptides_norm",
      values_fn = NULL,
    ) %>%
    tibble::column_to_rownames("quantCols")
  values <- heatmap_data %>%
    dplyr::select(-"condition", -"time") %>%
    as.matrix()
  colors <- heatmap_data %>%
    dplyr::select("condition", "time")
  peptide_order <- colnames(values)
  mat <- mspms_tidy_data %>%
    dplyr::select("quantCols", "peptide", "cleavage_seq") %>%
    tidyr::pivot_wider(
      names_from = "peptide",
      values_from = "cleavage_seq"
    ) %>%
    tibble::column_to_rownames("quantCols") %>%
    dplyr::select(dplyr::all_of(peptide_order)) %>%
    as.matrix()
  # Creating column labels for whether a peptide is cleaved or not
  col_logic <- grepl(".*_.*", colnames(values))
  col_p_colors <- c("cleavage_product", "full_length")
  col_colors <- data.frame(peptide_type = ifelse(col_logic,
    col_p_colors[1],
    col_p_colors[2]
  ))
  heatmaply::heatmaply(values,
    scale = scale,
    showticklabels = c(FALSE, TRUE),
    custom_hovertext = mat,
    row_side_colors = colors,
    col_side_colors = col_colors,
    plot_method = plot_method,
    show_dendrogram = show_dendrogram
  )
}
