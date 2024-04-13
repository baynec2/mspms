add_cleavages = function(joined_with_library){

  # Iterating through and applying nterm_clevage
  nterm = purrr::pmap_df(list(joined_with_library$Peptide,
                              joined_with_library$library_match_sequence,
                              joined_with_library$library_real_sequence),
                         mspms::nterm_cleavage)



  # Iterating though and applying cterm_cleavage
  cterm = purrr::pmap_df(list(joined_with_library$Peptide,
                              joined_with_library$library_match_sequence,
                              joined_with_library$library_real_sequence),
                         mspms::cterm_cleavage)



  # Combining nterm and cterm
  cleavages = dplyr::bind_cols(nterm,cterm[,2:3])


  # Building final data frame.
  output = dplyr::bind_cols(joined_with_library,cleavages) %>%
    dplyr::select(library_reference_id,library_real_sequence,Peptide,nterm,nterm_cleavage_pos,cterm,cterm_cleavage_pos,dplyr::everything(),-peptide)



  return(output)
}
