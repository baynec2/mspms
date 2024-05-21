#' plot_all_icelogos
#'
#' wrapper function for plotting icelogos for each condition in an experiment
#'
#' @param polished_data = this is the data that has been polished (via mspms::polish()) and is ready for plotting
#' @param type = this is the type of visualization you would like to perform, can be either "percent_difference" or "fold_change".
#' @param background_universe = this is a list of all the possible cleavage sequences in the peptide library used.
#'
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#'
#'
#' plot_all_icelogos(mspms::prepared_for_stats %>% mspms::polish())
plot_all_icelogos = function(polished_data,
                             type = "percent_difference",
                             background_universe = mspms::all_possible_8mers_from_228_library){

  # Extracting cleavage sequences
  cleavage_seqs = polished_data %>%
    dplyr::select(Peptide,cleavage_seq) %>%
    dplyr::distinct()


  # Generating stats, looking for an effect of time
  stats = mspms::mspms_anova(polished_data) %>%
    #extracting significant peptides
    dplyr::filter(p.adj <= 0.05) %>%
    dplyr::inner_join(cleavage_seqs,by = "Peptide")

  plot_list = list()

  # Plotting Icelogos
  for(i in unique(stats$condition)){
    #filtering data
    f = dplyr::filter(stats, condition == i) %>%
    #appending peptide to cleavage_seqs
      dplyr::pull(.data$cleavage_seq)

    # generating ice logo

  out = mspms::plot_icelogo(f,type = type,background_universe = background_universe)+
    ggplot::ggtitle(i)

  plot_list[[i]] = out

  }

  p1 =  do.call(ggpubr::ggarrange, c(plot_list, nrow = length(unique(stats$condition)),common.legend =TRUE,legend = "bottom"))

  return(p1)
}

