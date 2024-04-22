#' join_with_library
#'
#' Joins MSP-MS peptide library with imputed data.
#'
#' @param peptide_library = this is the peptide library that was used with the MSP-MS method
#' @param imputed = this is the imputed data. Intended to be prepared using mspms::impute()
#'
#' @return a data frame with the peptide library information joined to the imputed data.
#' @export
#'
#' @examples
#'
#' joined_with_library = join_with_library(mspms::imputed,mspms::peptide_library)
#'
join_with_library = function(imputed,peptide_library = mspms::peptide_library){

  # Joining the peptide library to the imputed data.
  joined = dplyr::inner_join(peptide_library,imputed,by = c("library_reference_id" = "Protein Accession")) %>%
    tibble::as.tibble()

  return(joined)
}
