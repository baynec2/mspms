#' calc_AA_count_of_motif
#'
#' calculate the counts of ammino acids at each position of a motif.
#'
#' @param cleavage_motif this is a vector of cleavage motifs that you would
#' like to count the amino acids for.
#' @return a matrix of counts
#' @export
#'
#' @examples
#' cleavage_motifs <- c("XXXXXXXX", "XXXXXXXX")
#' calc_AA_count_of_motif(cleavage_motifs)
#' #
calc_AA_count_of_motif <- function(cleavage_motif) {
  # calculating number of characters
  nchar <- nchar(cleavage_motif[1])

  # do this for the background first.
  sequences <- tibble::tibble(cleavage_motif = cleavage_motif)

  # separate the sequences into individual amino acids
  AAs <- tidyr::separate(
    sequences,
    col = 1,
    into = c("exclude", paste0("position", seq_len(nchar))),
    sep = "",
    remove = TRUE
  ) %>%
    dplyr::select(-.data$exclude)

  # Define the desired order.
  # This ensures that the amino acids will always be in the same order.
  desired_order <- c(
    "A", "D", "E", "F", "G", "H", "I", "K", "L", "n",
    "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y"
  )

  # Counting the number in each position
  count <- AAs %>%
    purrr::map_df(table)

  # Are any columns missing?
  `%!in%` <- Negate(`%in%`)
  missing_cols <- desired_order[desired_order %!in% names(count)]

  # building a df with missing columns
  missing_matrix <- matrix(nrow = nchar, ncol = length(missing_cols))
  colnames(missing_matrix) <- missing_cols
  missing_df <- tibble::as_tibble(missing_matrix)

  count_matrix <- dplyr::bind_cols(count, missing_df) %>%
    dplyr::relocate(desired_order) %>%
    t()

  # replacing NA with 0
  count_matrix[is.na(count_matrix)] <- 0

  colnames(count_matrix) <- paste0("P", c(
    seq(from = (nchar / 2), to = 1),
    paste0(seq_len(nchar / 2), "'")
  ))

  return(count_matrix)
}
