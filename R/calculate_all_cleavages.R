#' calculate_all_cleavages
#' calculate all possible cleavages for a defined peptide library containing
#' peptides of the same length.
#' @param peptide_library_seqs The sequences of each peptide in the peptide
#' library. They should all be the same length.
#' @param n_AA_after_cleavage  The number of AA after (and before) the cleavage
#' site to consider.
#' @return a vector of all the possible cleavages for the peptide library
#' sequences
#' @export
#' @examples
#' calculate_all_cleavages(mspms::peptide_library$library_real_sequence,
#'     n_AA_after_cleavage = 4
#' )
calculate_all_cleavages <- function(peptide_library_seqs,
                                    n_AA_after_cleavage = 4) {
    # calculate the length of the peptides in library.
    length <- unique(nchar(peptide_library_seqs))
    # The current code assumes all the sequences are the same length
    if (length(length) > 1) {
        stop("Peptide library sequences are not of the same length")
    }
    # first cut at every single position
    peptide_sequences <- c()
    # for each position go through and cut into two pieces
    for (i in seq_len(length - 1)) {
        ## Doing this for the first fragment, left side of cleavage
        first_fragment <- substr(peptide_library_seqs, 1, i)
        n_char <- unique(nchar(first_fragment))
        if (n_AA_after_cleavage - n_char > 0) {
            n_x <- paste0(rep("X", n_AA_after_cleavage - n_char), collapse = "")
        } else {
            n_x <- ""
        }
        first_fragment_x <- paste0(n_x, first_fragment)
        first_fragment_x <- substr(first_fragment_x,
            stop = unique(nchar(first_fragment_x)),
            start = unique(nchar(first_fragment_x)) -
                n_AA_after_cleavage + 1
        )
        ## Doing this for the second fragment, right side of cleavage
        second_fragment <- substr(
            peptide_library_seqs, i + 1,
            i + n_AA_after_cleavage
        )
        n_char <- unique(nchar(second_fragment))
        # calculating the overhangs
        if (n_AA_after_cleavage - n_char > 0) {
            n_x <- paste0(rep("X", n_AA_after_cleavage - n_char), collapse = "")
        } else {
            n_x <- ""
        }
        # appending the overhangs
        second_fragment_x <- paste0(second_fragment, n_x)
        # combining the two fragments
        combined <- paste0(first_fragment_x, second_fragment_x)
        # appending the combined fragments to the vector
        peptide_sequences <- c(peptide_sequences, combined)
    }
    return(peptide_sequences)
}
