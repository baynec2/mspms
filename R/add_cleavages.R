#' add_cleavages
#'
#' This functions adds cleavage information to the tibble that has been normalized, outlier removed, imputed, and joined with (peptide) library.
#' The cleavage information is represented as the 4 amino acids to the left and right of the cleavage site.
#' If there is no amino acid in a position of original peptide in the library that was cleaved, it is represented as an X.
#' This wraps the mspms::n_term_cleavage and mspms::c_term_cleavage functions into a consolidated function.
#'
#' @param joined_with_library = this is the tibble that has been normalized, outlier removed, imputed, and joined with the library.
#' @param n_residues = the number of residues to the left and right of the cleavage site to include in the output.
#' @return a tibble with cleavage information added.
#' @export
#'
#' @examples
#' # adding the clevages 5 AA to the left and right of the cleavage site.
#' add_cleavages(mspms::joined_with_library,n_residues = 5)
#'
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
    dplyr::select(.data$Peptide,
                  .data$library_reference_id,
                  .data$library_match_sequence,
                  .data$library_real_sequence,
                  .data$nterm,
                  .data$nterm_cleavage_pos,
                  .data$cterm,
                  .data$cterm_cleavage_pos,
                  dplyr::everything(),
                  -.data$peptide,
                  -.data$Peptide_no_cleavage,
                  -dplyr::any_of("RT")) %>%
    tibble::as_tibble()



  return(output)
}
