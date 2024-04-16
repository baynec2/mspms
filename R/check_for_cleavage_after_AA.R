#' check_for_cleavage_after_AA
#'
#' This function checks if the cleavage site is after a vector of specified amino acids.
#'
#' Note that the cleavage sequence is intended to be represented as 4 AA to the left and right of cleavage site, as single letter codes.
#'
#' @param cleavage_seq = this is the cleavage sequecne that is being checked. Represented as 4 AA to the left and right of cleavage site
#' @param AA_vector = this is a vector of AA that you want to check to see if the cleavage site is after. For example c("W","L","Y","F").
#'
#' @return
#' a logical vector of the same length as the input cleavage_seq
#' @export
#'
#' @examples
#'
#' #this is true because the cleavage site is after a W
#' check_for_cleavage_after_AA("GLYWAFAA",c("W","L","Y","F"))
#'
#' #this is false because the cleavage site is not after a W, L, Y, or F
#' check_for_cleavage_after_AA("GLYAAFAA",c("W","L","Y","F"))
#'
check_for_cleavage_after_AA = function(cleavage_seq,AA_vector){

  # check if the cleavage site is in the AA_list
  out = stringr::str_sub(cleavage_seq,4,4) %in% AA_list

  return(out)
}
