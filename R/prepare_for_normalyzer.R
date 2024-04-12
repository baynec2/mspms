#' prepare_for_normalyzer: This function prepares the data out from peaks for downstream normalization via normalyzer.
#' 
#' Briefly, it combines the two files by peptide IDs and filters out peptides that have quality scores <= 0.3. 0s are also replaced with NAs 
#' Note that there can be multiple matches between these, currently  only keeps the values for the best match. This could be a bug. 
#' @param lfq_filepath = this is the filepath to the first PEAKS output table
#' @param id_filepath = this is the filepath to the second PEAKS output table
#'
#' @return a data frame containing the combined columns from the 
#' @export
#'
#' @examples
prepare_for_normalyzer = function(lfq_filepath,
                                  id_filepath){
  
  # Reading in the label free quantification data  
  lfq = readr::read_csv(lfq_filepath,guess_max = 10)
  
  # Reading in the ids
  id = readr::read_csv(id_filepath) %>% 
    #only keep the Peptide ID with the highest score in case there are more than one. 
    # This is essentially what the current script does, is it intended? Seems like a bug to me
    dplyr::group_by(Peptide) %>% 
    dplyr::filter(`-10lgP` == max(`-10lgP`)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(Peptide,RT,6:12)

  
  # Combining data frame and filtering to only contain quality scores > 0.3
  output = lfq %>% 
    dplyr::inner_join(id,by = c("Peptide"),multiple = "first") %>% 
    dplyr::filter(Quality > 0.3)
  
  # Replacing 0 with NA
  output[output == 0] = NA_real_
  
  return(output)
}
