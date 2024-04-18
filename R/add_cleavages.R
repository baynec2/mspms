#' add_cleavages
#'
#' This functions adds cleavage information to the data frame that has been normalized, outlier removed, imputed, and joined with (peptide) library.
#' The cleavage information is represented as the 4 amino acids to the left and right of the cleavage site.
#' If there is no amino acid in a position of original peptide in the library that was cleaved, it is represented as an X.
#' This wraps the mspms::n_term_cleavage and mspms::c_term_cleavage functions into a consolidated function.
#'
#' @param joined_with_library
#'
#' @return
#' a data frame with cleavage information added.
#' @export
#'
#' @examples
add_cleavages = function(joined_with_library,n_residues = 4){

  # Iterating through and applying nterm_clevage
  nterm = purrr::pmap_df(list(joined_with_library$Peptide,
                              joined_with_library$library_match_sequence,
                              joined_with_library$library_real_sequence,
                              n_residues),
                         mspms::nterm_cleavage)



  # Iterating though and applying cterm_cleavage
  cterm = purrr::pmap_df(list(joined_with_library$Peptide,
                              joined_with_library$library_match_sequence,
                              joined_with_library$library_real_sequence,
                              n_residues),
                         mspms::cterm_cleavage)



  # Combining nterm and cterm
  cleavages = dplyr::bind_cols(nterm,cterm[,2:3])


  # Building final data frame.
  output = dplyr::bind_cols(joined_with_library,cleavages) %>%
    dplyr::select(library_reference_id,library_real_sequence,Peptide,nterm,nterm_cleavage_pos,cterm,cterm_cleavage_pos,dplyr::everything(),-peptide) %>%
    tibble::as_tibble()



  return(output)
}
