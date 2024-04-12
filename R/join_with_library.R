#' join_with_library
#'
#' @param peptide_library
#' @param imputed
#'
#' @return
#' @export
#'
#' @examples
join_with_library = function(imputed,peptide_library = mspms::peptide_library){

  # Joining the peptide library to the imputed data.
  joined = dplyr::inner_join(peptide_library,imputed,by = c("library_reference_id" = "Protein Accession"))

  return(joined)
}
