#' plot_icelogo
#'
#' This function plots the cleavage motifs that were enriched relative to background in a manner similar to IceLogo.
#' https://iomics.ugent.be/icelogoserver/resources/manual.pdf
#'
#' @param cleavage_seqs = these are the cleavage sequences that are observed in the experiment
#' @param background_universe = this is a list of all the possible cleavage sequences in the peptide library used.
#' @param pval = this is the pvalue threshold that you would like to use.
#' @param type = this is the type of visualization you would like to perform, can be either "percent_difference" or "fold_change".
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#' cleavage_seqs <- mspms::mspms_data %>%
#'   dplyr::filter(condition == "DMSO", time == 240) %>%
#'   dplyr::pull(cleavage_seq)
#'
#' plot_icelogo(cleavage_seqs)
plot_icelogo <- function(cleavage_seqs,
                         background_universe = mspms::all_possible_8mers_from_228_library,
                         pval = 0.05,
                         type = "percent_difference") {
  # calculation proportions of background
  background_counts <- mspms::calc_AA_count_of_motif(background_universe)
  background_proportions <- mspms::calc_AA_prop_of_motif(background_counts)

  # calculate proportions of experimentally identified motifs
  experimental_counts <- mspms::calc_AA_count_of_motif(cleavage_seqs)
  experimental_proprotions <- mspms::calc_AA_prop_of_motif(experimental_counts)

  # calculate zscores
  zscores <- mspms::calc_AA_motif_zscore(
    background_count_matrix = background_counts,
    background_prop_matrix = background_proportions,
    experimental_count_matrix = experimental_counts,
    experimental_prop_matrix = experimental_proprotions
  )

  # calculate significant zscores
  sig_zscores <- mspms::calc_sig_zscores(zscores, pval)

  if (type == "percent_difference") {
    pd <- mspms::calc_AA_percent_difference(
      experimental_prop_matrix = experimental_proprotions,
      background_prop_matrix = background_proportions
    )
    final_pd <- mspms::prepare_sig_p_dif(pd, sig_zscores)
    plot <- mspms::plot_pd_icelogo(final_pd)
    return(plot)
  } else if (type == "fold change") {
    fc <- mspms::calc_AA_fc(
      experimental_prop_matrix = experimental_proprotions,
      background_prop_matrix = background_proportions
    )
    final_fc <- mspms::prepare_fc(fc, sig_zscores)
    plot <- mspms::plot_fc_icelogo(final_fc)
    return(plot)
  } else {
    stop("type must be either 'percent_difference' or 'fold change'")
  }
}
