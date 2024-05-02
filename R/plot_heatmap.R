#' plot_heatmap
#'
#' This produces a heatmaply interactive heatmap of the data with color bars representing the condition and time for each sample.
#'
#' @param prepared_for_stats = this is the data that has been processed by the prepare for stats function.
#' @param scale = how would you like the data scaled? default is none, but can be done by "row", "column", or "none"
#'
#' @return a heatmaply interactive heatmap
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' plot_heatmap(mspms::prepared_for_stats,scale = "column")
#'
plot_heatmap = function(prepared_for_stats,scale = "column"){

  heatmap_data = prepared_for_stats %>%
    dplyr::select(.data$sample,.data$Peptide,.data$condition,.data$time,.data$value) %>%
    tidyr::pivot_wider(names_from = .data$Peptide,values_from = .data$value,values_fn = mean) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::mutate(time = as.factor(.data$time))


 heatmaply::heatmaply(heatmap_data,scale = scale,showticklabels = c(FALSE,TRUE))


}


