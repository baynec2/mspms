#' plot_heatmap
#'
#' This produces a heatmaply interactive heatmap of the data with color bars representing the condition and time for each sample.
#'
#' @param mspms_data = this is the data that has been processed by the mspms pipeline.
#' @param scale = how would you like the data scaled? default is none, but can be done by "row", "column", or "none"
#'
#' @return a heatmaply interactive heatmap
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' plot_heatmap(mspms::prepared_for_stats,scale = "column")
#'
plot_heatmap = function(mspms_data,scale = "column"){


  heatmap_data = mspms_data %>%
    dplyr::select(sample,Peptide,condition,time,value) %>%
    tidyr::pivot_wider(names_from = Peptide,values_from = value,values_fn = mean) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::mutate(time = as.factor(time))


  peptide_order = names(heatmap_data)[3:length(heatmap_data)]

  lab = mspms_data %>%
    dplyr::select(sample,Peptide,cleavage_seq,condition,time) %>%
    tidyr::pivot_wider(names_from = Peptide,values_from = cleavage_seq) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(condition,time,dplyr::all_of(peptide_order)) %>%
    as.matrix()

 heatmaply::heatmaply(heatmap_data,
                      scale = scale,
                      showticklabels = c(FALSE,TRUE),
                      custom_hovertext = lab
              )


}



