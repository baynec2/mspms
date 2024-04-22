#' plot_cleavage_motifs
#'
#' This function plots the cleavage motifs that were enriched relative to background in a manner similar to IceLogo.
#' https://iomics.ugent.be/icelogoserver/resources/manual.pdf
#'
#' @param cleavage_seqs = these are the cleavage sequences that are observed in the experiment
#' @param background_universe = this is a list of all the possible clevage sequences in the peptide library used.
#'
#' @return a ggplot object that shows the motif of the cleavage sequences
#' @export
#'
#' @examples
#'
#'cleavage_seqs = mspms::cleavage_added_data %>%
#'filter(condtion == "DMSO",time == 240) %>%
#'pull(cleavage_seq)
#'
#'background_universe = mspms::mspms::all_possible_8mers_from_228_library
#'
#'
#'plot_clavage_motif(cleavage_seqs,background_universe)

plot_cleavage_motif = function(cleavage_seqs,background_universe = mspms::all_possible_8mers_from_228_library){

   ## First need to format this such that we have a matrix with the counts at each location of the sequence
  nchar_bg = nchar(background_universe[1])

  # do this for the background first.
  background_universe = tibble(sequences = background_universe)

  # separate the sequences into individual amino acids
  bg_seq = tidyr::separate(background_universe,col = 1,
                           into = paste0("position",1:(nchar_bg+1)),
                           sep = "",remove = TRUE) %>%
    dplyr::select(-position1)

  # Counting the number in each position
  bg_count = bg_seq %>%
    purrr::map_df(table) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("AA")

  #Fixing names
  names(bg_count) = c("AA",paste0("position",1:nchar_bg))

  bg_prop = bg_count %>%
    mutate(across(where(is.numeric), prop.table)) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()

  bg_bits = mspms::calc_bits(bg_prop)

 # Now doing this for the experimental clevages

 nchar_cleav = nchar(cleavage_seqs[[1]])

 # done for the background first
 cleavage_seqs = tibble(sequences = cleavage_seqs)

 clev_seq = tidyr::separate(cleavage_seqs,col = 1,
                            into = paste0("position",1:(nchar_cleav+1)),
                            sep = "",remove = TRUE) %>%
   dplyr::select(-position1)

 # counting the number of time each AA appears at each position

 clev_count = clev_seq %>%
   purrr:::map_df(table) %>%
   t() %>%
   as.data.frame() %>%
   tibble::rownames_to_column("AA")

  clev_count[is.na(clev_count)] = 0
  names(clev_count) = c("AA",paste0("position",1:nchar_cleav))

  clev_prop = clev_count %>%
   mutate(across(where(is.numeric), prop.table)) %>%
   tibble::column_to_rownames("AA") %>%
   as.matrix()

  # calculating percentage difference between experimenally observed and theoritcal cleavage propensities
  clev_bits = clev_prop %>%
    mspms::calc_bits()

  # Subtracting the background bits from the cleavage bits
  bg_sub_bits = clev_bits - bg_bits

# Plotting the motif
p1 = ggseqlogo::ggseqlogo(bg_sub_bits,font = "helvetica_light" , method='custom',seq_type = "AA") +
  ggplot2::scale_x_continuous(labels = paste0("P",c((nchar_cleav/2):1,
                                                    paste0(1:(nchar_cleav/2),"'"))),
                              breaks = 1:nchar_cleav)
return(p1)
}



