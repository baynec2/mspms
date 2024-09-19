#' plot_heatmap
#'
#' This produces a heatmaply interactive heatmap of the QFeatures object with 
#' color bars representing the condition and time for each sample in each row.
#' 
#' Each column has a colored bar representing whether the peptide is a 
#' cleavage product or a full length member of the peptide library.
#'
#' @param processed_qf a QFeatures object containing a SummarizedExperiment
#' named "peptides_norm"
#' @param scale how would you like the data scaled? default is none,
#'  but can also be "row", "column", or "none"
#' @param plot_method what plot method would you like to use, can use
#' plotly or ggplot2. 
#'
#' @return a heatmaply interactive heatmap
#' @export
#' 
#' @examples
#' plot_heatmap(mspms::processed_qf)

plot_heatmap <- function(processed_qf,
                         scale = "column",
                         plot_method = "plotly") {
  # tidying the data
  mspms_data <- mspms_tidy(processed_qf,"peptides_norm")
  heatmap_data <- mspms_data %>%
    dplyr::select("quantCols","condition","time","peptide", "peptides_norm") %>%
    tidyr::pivot_wider(
      names_from = "peptide",
      values_from = "peptides_norm",
      values_fn = NULL,
    ) %>%
    tibble::column_to_rownames("quantCols") 
    
    values = heatmap_data %>% 
      dplyr::select(-"condition",-"time") %>% 
      as.matrix() 
    
    colors = heatmap_data %>% 
      dplyr::select("condition","time")
    
    peptide_order <- colnames(values)
    
    mat <- mspms_data %>%
      dplyr::select("quantCols", "peptide", "cleavage_seq") %>%
      tidyr::pivot_wider(names_from = "peptide",
                         values_from = "cleavage_seq") %>%
      tibble::column_to_rownames("quantCols") %>%
      dplyr::select(dplyr::all_of(peptide_order)) %>%
      as.matrix()
    # Creating column labels for whether a peptide is cleaved or not
    col_logic = grepl(".*_.*",colnames(values))
    col_p_colors = c("cleavage_product","full_length")
    col_colors = data.frame(peptide_type = ifelse(col_logic,
                                      col_p_colors[1],
                                      col_p_colors[2]))
    heatmaply::heatmaply(values,
      scale = scale,
      showticklabels = c(FALSE, TRUE),
      custom_hovertext = mat,
      row_side_colors = colors,
      col_side_colors = col_colors,
      plot_method = plot_method)
  
}
