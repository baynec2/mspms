#' Plot cleavage-site amino acid frequencies for an experimental set relative
#' to a reference set
#'
#' This function calculates amino acid frequencies at each motif position for
#' an experimental cleavage sequence set, compares them to a reference set,
#' filters residues based on a p-value threshold, and returns an icelogo-style
#' bar plot.
#'
#' @param cleavage_seqs Character vector of experimental cleavage sequences.
#' @param background_universe Character vector representing the reference
#'   cleavage sequence set. Defaults to all possible 8-mers from the 228
#'   library.
#' @param pval_threshold Numeric value. Residues with p-values greater than
#'   this threshold will be omitted from the plot.
#'
#' @return A ggplot2 object showing amino acid frequencies at each motif
#'   position for the experimental and reference sets.
#' @export
#'
#' @examples
#' reference_set = rep(c("AAAAAAAA"),10)
#' plot_cleavage_frequencies(reference_set)
plot_cleavage_frequencies <- function(cleavage_seqs,
                                      background_universe =
                                        mspms::all_possible_8mers_from_228_library,
                                      pval_threshold = 1) {
  if (length(cleavage_seqs) != length(unique(cleavage_seqs))) {
    message("non unique cleavage sequences supplied. Only considering unique
            sequences")
    cleavage_seqs <- unique(cleavage_seqs)
  }

  # Frequencies of each
  background_counts <- mspms:::calc_AA_count_of_motif(background_universe)
  background_proportions <- mspms:::calc_AA_prop_of_motif(background_counts)
  experimental_counts <- mspms:::calc_AA_count_of_motif(cleavage_seqs)
  experimental_proprotions <- mspms:::calc_AA_prop_of_motif(experimental_counts)

  # Pvalues
  zscores <- mspms:::calc_AA_motif_zscore(
    background_count_matrix = background_counts,
    background_prop_matrix = background_proportions,
    experimental_count_matrix = experimental_counts,
    experimental_prop_matrix = experimental_proprotions
  )
  # Converting Zscores into pvalues
  pvalues <- 2 * (1 - pnorm(abs(as.matrix(zscores))))

  # Assemble into dataframe
  background <- as.data.frame(background_proportions) %>%
    tibble::rownames_to_column("amino_acid") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "frequency") %>%
    dplyr::mutate(set = "reference")

  experimental <- as.data.frame(experimental_proprotions) %>%
    tibble::rownames_to_column("amino_acid") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "frequency") %>%
    dplyr::mutate(set = "experimental")

  pvalues <- as.data.frame(pvalues) %>%
    tibble::rownames_to_column("amino_acid") %>%
    tidyr::pivot_longer(2:length(.), names_to = "position", values_to = "pvalues")

  # Combine Everything
  data <- dplyr::bind_rows(background, experimental)
  data_pval <- dplyr::inner_join(data, pvalues, by = c("amino_acid", "position"))

  # Fixing levels/
  levels(data_pval$position) <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")

  p1 <- data_pval %>%
    ggplot2::ggplot(ggplot2::aes(amino_acid, frequency, fill = set)) +
    ggplot2::geom_col(position = "dodge")

  if (pval_threshold < 1) {
    sig_annotation <- ggplot2::geom_text(
      data <- data_pval %>%
        dplyr::group_by(position, amino_acid) %>%
        dplyr::filter(frequency == max(frequency)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(pvalues <= pval_threshold),
      mapping = ggplot2::aes(
        x = amino_acid,
        y = frequency,
        label = "*"
      ),
      color = "black",
      size = 6,
      inherit.aes = FALSE
    )
    p1 <- p1 + sig_annotation
  }
  # Getting the facets right
  p1 <- p1 +
    ggplot2::facet_wrap(~position, ncol = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::theme(legend.position = "bottom")

  return(p1)
}
