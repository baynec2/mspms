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
  # calculate n 
  n = length(cleavage_seqs)
  
# processing data to plot
 icelogo_data = mspms::prepare_icelogo_data(cleavage_seqs,
                               background_universe,
                               pval,
                               type)
  
 nchar_cleav <- ncol(icelogo_data)
 
 # Plotting the motif
 p1 <- ggseqlogo::ggseqlogo(icelogo_data,
                            font = "helvetica_light",
                            method = "custom",
                            seq_type = "AA") +
   ggplot2::scale_x_continuous(
     labels = paste0("P", c(
       (nchar_cleav / 2):1,
       paste0(
         seq_len(nchar_cleav / 2), "'"
       )
     )),
     breaks = seq_len(nchar_cleav)
   ) +
   ggplot2::labs(subtitle = paste0("n = ",n))+
   if(type == "percent_difference") {
     ggplot2::ylab("percent difference")
   } else if(type == "fold_change") {
     ggplot2::ylab("fold change")
   }else{
     stop("type must be either percent_difference or fold_change")
   }
 
 return(p1)
}

