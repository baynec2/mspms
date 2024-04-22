#' calc_bits
#'
#'function to calculate the bits of a given proportion as described: https://iomics.ugent.be/icelogoserver/resources/manual.pdf
#'
#' @param proportion  a numeric value representing the proportion of a given amino acid
#'
#' @return a numeric value representing the bits of the given proportion
#' @export
#'
#' @examples
#' calc_bits(0.25)
calc_bits = function(proportion){
  bits = (proportion * log2(proportion)) * -1
  return(bits)
}
