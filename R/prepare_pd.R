#' prepare_pd
#'
##' Prepare a label free quantification file exported from Proteome Discoverer
#' for subsequent mspms analysis.
#'
#' @param filepath = this is the filepath to the proteome discoverer excel
#' formatted file.
#'
#' @return a tibble with the data formatted for use with normalyze
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#' prepared_proteome_discoverer <- prepare_pd("tests/testdata/proteome_discoverer_output.xlsx")
prepare_pd <- function(filepath) {
  
  return(sample_name_fixed)
}
