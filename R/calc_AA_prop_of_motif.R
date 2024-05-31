#' calc_AA_prop_of_motif
#'
#' calculate the proportion of amino acids at each position in a vector of
#' motifs.
#'
#' @param count_matrix this is a matrix of the counts of cleavage motifs
#'
#' @return a matrix with proportions of counts.
#' @export
#' @examples
#'
#' count_matrix <- matrix(c(1, 10, 100, 1000, 9, 90, 900, 9000),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#' calc_AA_prop_of_motif(count_matrix)
#'
calc_AA_prop_of_motif <- function(count_matrix) {
  prop <- count_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("AA") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), prop.table)) %>%
    tibble::column_to_rownames("AA") %>%
    as.matrix()


  return(prop)
}
