#' calc_AA_motif_zscore
#'
#' calculate the z score for AAs at each position
#'
#' @param background_count_matrix = this is the count matrix from the
#' background sequences
#' @param background_prop_matrix = this is the proportion matrix from the
#' background sequences
#' @param experimental_count_matrix = this is the count matrix from the
#' experimental sequences
#' @param experimental_prop_matrix = this is the proportion matrix from the
#' experimental sequences
#'
#' @return a data frame of zscores for each AA at each position.
#' @export
#' @examples
#' background_count_matrix <- matrix(c(10, 10, 10, 10),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' background_prop_matrix <- matrix(c(0.4, 0.4, 0.4, 0.4),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' experimental_count_matrix <- matrix(c(5, 5, 5, 5),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' experimental_prop_matrix <- matrix(c(0.2, 0.2, 0.2, 0.2),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' calc_AA_motif_zscore(
#'   background_count_matrix,
#'   background_prop_matrix,
#'   experimental_count_matrix,
#'   experimental_prop_matrix
#' )
calc_AA_motif_zscore <- function(background_count_matrix,
                                 background_prop_matrix,
                                 experimental_count_matrix,
                                 experimental_prop_matrix) {
  ### calculating the SD ###
  # defining function per Icelogo manual
  standard_deviation <- function(proportion_matrix, count_matrix) {
    (proportion_matrix / colSums(count_matrix))^(1 / 2)
  }

  # calculating the standard deviation
  bg_sd <- standard_deviation(
    background_prop_matrix,
    experimental_count_matrix
  )

  ### Calculating Z score ###

  # Defining function per icelogo manual
  zscore <- function(experimental_prop_matrix, background_prop_matrix, bg_sd) {
    (experimental_prop_matrix - background_prop_matrix) / bg_sd
  }

  # applying to our data
  zscores <- zscore(experimental_prop_matrix, background_prop_matrix, bg_sd) %>%
    as.data.frame()

  return(zscores)
}
