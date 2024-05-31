#' calc_AA_percent_difference
#'
#' calculate the percent difference between a matrix of background proportions
#' and a matrix of experimentally observed proportions.
#'
#' @param background_prop_matrix = proportion matrix of aa per position from
#'  background cleavage sequences
#' @param experimental_prop_matrix = a proportion matrix of aa per position from
#'  experimental cleavage sequences
#'
#' @return a data frame of percent differences
#' @export
#' @examples
#' experimental_prop_matrix <- matrix(c(0.2, 0.2, 0.2, 0.2),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' background_prop_matrix <- matrix(c(0.4, 0.4, 0.4, 0.4),
#'   nrow = 1,
#'   dimnames = list("A")
#' )
#' calc_AA_percent_difference(background_prop_matrix, experimental_prop_matrix)
calc_AA_percent_difference <- function(background_prop_matrix,
                                       experimental_prop_matrix) {
  # Calculating percent difference
  pd <- (experimental_prop_matrix * 100) - (background_prop_matrix * 100) %>%
    as.data.frame()

  return(pd)
}
