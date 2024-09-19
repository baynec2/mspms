#' plot_all_icelogos
#'
#' Easily plot a iceLogo corresponding to peptides of interest across each 
#' condition of an experiment.
#' @param sig_cleavage_data a tibble of data of interest containing a column
#' labeled peptide, cleavage_seq, and condition
#' @param pval this is the pvalue threshold (<=) to consider significant when 
#' determining the significance of the sig_cleavages relative to the background 
#' at each position of the iceLogo.
#' @param typethis is the type of iceLogo you would like to generate,
#'  can be either "percent_difference" or "fold_change".
#' @param background_universe this is a list cleavages you would like to compare
#' to as background of the iceLogo
#'
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#' # Determining cleavages of interest
#' sig_cleavage_data <- mspms::log2fc_t_test_data %>% 
#' dplyr::filter(p.adj <= 0.05,log2fc > 3) 
#' # Plotting a iceLogo for each condition.
#'plot_all_icelogos(sig_cleavage_data)
plot_all_icelogos <- function(sig_cleavage_data,
                              type = "percent_difference",
                              pval = 0.05,
                              background_universe =
                                mspms::all_possible_8mers_from_228_library) {
  # initiating plot list
  plot_list <- list()

  # Plotting Icelogos
  for (i in unique(sig_cleavage_data$condition)) {
    # filtering data
    f <- dplyr::filter(sig_cleavage_data, .data$condition == i) %>%
      dplyr::select("peptide", "cleavage_seq") %>%
      # appending peptide to cleavage_seqs
      dplyr::pull(.data$cleavage_seq) %>% 
      unique()

    # generating ice logo
    out <- mspms::plot_icelogo(f,
      type = type, pval = pval,
      background_universe = background_universe
    ) +
      ggplot2::ggtitle(i)

    plot_list[[i]] <- out
  }
  p1 <- do.call(ggpubr::ggarrange, c(plot_list,
    nrow = length(unique(sig_cleavage_data$condition)),
    common.legend = TRUE, legend = "bottom"
  ))
  return(p1)
}
