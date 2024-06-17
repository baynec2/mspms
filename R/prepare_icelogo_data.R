#' prepare_icelogo_data
#' 
#' prepare a matrix containing icelogo data. 
#' 
#' @param cleavage_seqs = these are the cleavage sequences that are observed in the experiment
#' @param background_universe = this is a list of all the possible cleavage sequences in the peptide library used.
#' @param pval  = this is the pvalue threshold that you would like to use.
#' @param type = this is the type of Icelogo matrix you would like to prepare, can be either "percent_difference" or "fold_change".
#'
#' @return a matrix of enriched amino acids per position
#' @export
#'
#' @examples
#' cleavage_seqs = unique(mspms::mspms_data$cleavage_seq[1:100])
#' prepare_icelogo_data(cleavage_seqs)
#' 
prepare_icelogo_data <- function(cleavage_seqs,
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
    return(final_pd)
  } else if (type == "fold_change") {
    fc <- mspms::calc_AA_fc(
      experimental_prop_matrix = experimental_proprotions,
      background_prop_matrix = background_proportions,
      sig_zscores = sig_zscores
    )
    return(fc)
  } else {
    stop("type must be either 'percent_difference' or 'fold_change'")
  }
}
