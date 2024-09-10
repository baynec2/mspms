#' plot_all_icelogos
#'
#' wrapper function for plotting icelogos for each condition in an experiment.
#'  Here we will plot all icelogos for all peptides that were found to be
#'  significant by t test at any timepoint.
#'
#' @param mspms_data = this is the data that has been polished
#' (via mspms::polish()) and is ready for plotting
#' @param pval = this is the pvalue threshold that you would like to use.
#' @param log2fc_threshold = this is the log2 fold change threshold that you 
#' would like to use. Can be positive or negative if desired. 
#' @param type = this is the type of visualization you would like to perform,
#'  can be either "percent_difference" or "fold_change".
#' @param background_universe = this is a list of all the possible cleavage
#' sequences in the peptide library used.
#'
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#' plot_all_icelogos(mspms::mspms_data)
plot_all_icelogos <- function(mspms_data,
                              pval = 0.05,
                              log2fc_threshold = 3,
                              type = "percent_difference",
                              background_universe =
                                mspms::all_possible_8mers_from_228_library) {
  # Generating stats,

  if(log2fc_threshold >= 0){
    stats <- mspms::log2fc_t_test(mspms_data) %>%
      # extracting significant peptides
      dplyr::filter(
        .data$p.adj <= pval,
        .data$log2fc > log2fc_threshold
      )
  } else {
    stats <- mspms::log2fc_t_test(mspms_data) %>%
      # extracting significant peptides
      dplyr::filter(
        .data$p.adj <= pval,
        .data$log2fc < log2fc_threshold & .data$log2fc < 0
      )
  }
  plot_list <- list()

  # Plotting Icelogos
  for (i in unique(stats$condition)) {
    # filtering data
    f <- dplyr::filter(stats, .data$condition == i) %>%
      dplyr::select("Peptide", "cleavage_seq") %>%
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
    nrow = length(unique(stats$condition)),
    common.legend = TRUE, legend = "bottom"
  ))

  return(p1)
}
