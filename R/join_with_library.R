#' join_with_library
#'
#' @param peptide_library 
#' @param imputed 
#'
#' @return
#' @export
#'
#' @examples
join_with_library = function(imputed,peptide_library = readRDS("peptide_library")){
  
  # Joining the peptide library to the imputed data. 
  joined = dplyr::right_join(peptide_library,imputed,by = c("library_match_sequence" = "Peptide_no_cleavage"))
  
  return(joined)
}
