#' plot_heatmap
#'
#' This produces a heatmaply interactive heatmap of the data with color bars
#'  representing the condition and time for each sample.
#'
#' @param mspms_data = this is the data that has been processed by the
#' mspms pipeline.
#' @param scale = how would you like the data scaled? default is none,
#'  but can be done by "row", "column", or "none"
#'
#' @return a heatmaply interactive heatmap
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' plot_heatmap(mspms::mspms_data, scale = "column")
#'
plot_heatmap <- function(mspms_data, scale = "column") {
  
  heatmap_data <- mspms_data %>%
    dplyr::select("sample","condition","time","Peptide", "value") %>%
    tidyr::pivot_wider(
      names_from = "Peptide",
      values_from = "value",
      values_fn = NULL,
    ) %>%
    tibble::column_to_rownames("sample") 
    
    values = heatmap_data %>% 
      dplyr::select(-"condition",-"time") %>% 
      as.matrix() 
    
    colors = heatmap_data %>% 
      dplyr::select("condition","time")
    
    peptide_order <- colnames(values)
    
    mat <- mspms_data %>%
      dplyr::select("sample", "Peptide", "cleavage_seq") %>%
      tidyr::pivot_wider(names_from = "Peptide",
                         values_from = "cleavage_seq") %>%
      tibble::column_to_rownames("sample") %>%
      dplyr::select(dplyr::all_of(peptide_order)) %>%
      as.matrix()
    # Creating column labels for whether a peptide is cleaved or not
    col_logic = grepl(".*_.*",colnames(values))
    col_p_colors = c("cleaved","not_cleaved")
    col_colors = data.frame(cleavage_status = ifelse(col_logic,
                                      col_p_colors[1],
                                      col_p_colors[2]))
    heatmaply::heatmaply(values,
      scale = scale,
      showticklabels = c(FALSE, TRUE),
      custom_hovertext = mat,
      row_side_colors = colors,
      col_side_colors = col_colors
    )
  
}
