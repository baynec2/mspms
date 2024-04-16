#' cterm_cleavage
#'
#' Finding the cleavage sequences on the C terminus of a given peptide in reference to the peptide library it was derived from
#'
#' @param peptide_sequence = this is the peptide sequence in single leter AA code. _ denotes cleavage site.
#' @param library_match_sequence  = this is the sequence of the corresponding peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the Ls transformed to Ms (This is what the legacy code did so it is kept this way in case there was a good reason for it)
#'
#' @return
#' a data frame with the peptide sequence, cleavage sequences 4 AA on the left and right of the c term cleavage, and the position of the c term cleavage in the library sequence.
#' @export
#'
#' @examples
#' peptide_sequence = "YWLSTHLAGKR_R"
#' library_match_sequence = "YWLSTHLAGKRRDW"
#' library_real_sequence = "YWMSTHLAGKRRDW"
#' cterm_cleavage(peptide_sequence,library_match_sequence,library_real_sequence)
#'
cterm_cleavage = function(peptide_sequence,library_match_sequence,library_real_sequence){

  # First, we determine if there is a cterm cleavage
  # We do this by looking to see if there is a _ (denotes a cleavage)
  # If there is, we check to see if it is in the second to last position of the peptide when there is only one _ (only cterm cleavage)
  # Or there might be 2 _s in the peptide sequence (meaning there is a cterm and nterm cleavage)

  n_of_match = length(gregexpr('_',peptide_sequence)[[1]])
  if((grepl("_", peptide_sequence)==TRUE) && (gregexpr('_',peptide_sequence)[[1]][n_of_match]==nchar(peptide_sequence)-1)){

    # if there is a c term cleavage, it is the last position - 1 of the peptide sequence
    pos = nchar(peptide_sequence) - 1

    # now we define the sequence on the left side of the cleavage
    temp = substr(peptide_sequence,pos-4, pos-1)

    # Checking to see what part of the reference sequence this matches.
    left_reference_beginning = regexpr(temp,library_match_sequence)[[1]][1]
    # it is four AA long, so the end is the beginning + 3
    left_reference_end = left_reference_beginning + 3

    # now we can determine the right side of the cleavage sequence by taking the left reference end and adding 1
    right_reference_beginning = left_reference_end + 1

    # the cterm cleavage position is the right reference beginning minus 1
    cterm_cleavage_pos = right_reference_beginning -1


    # it is four AA long, so the end is the beginning + 3
    right_reference_end = right_reference_beginning + 3

    # Deciding what to do based on what position of the reference the match is to.
    # Note that all of the peptides in the library are 14 AA long
    # This is all about getting the 4 flanking amino acids correct. It should be 4 on each side of the cleavage
    # X denotes that there was nothing there, as happens if it is on the edge of the sequence in the peptide library.
    if(left_reference_end == 13){
      cterm = paste(
        c(substr(library_real_sequence,left_reference_beginning,left_reference_end),
          substr(library_real_sequence,right_reference_beginning,right_reference_beginning)
          ,"XXX"),collapse = "")
    }else if(left_reference_end == 12){
      cterm = paste(
        c(
          substr(library_real_sequence,left_reference_beginning,left_reference_end),
          substr(library_real_sequence,right_reference_beginning,right_reference_beginning+1),
                 "XX"),collapse = "")
    }else if(left_reference_end == 11){
      cterm = paste(
        c(
          substr(library_real_sequence,left_reference_beginning,left_reference_end),
          substr(library_real_sequence,right_reference_beginning,right_reference_beginning +2),
          "X"),collapse = "")
    }else if(dplyr::between(left_reference_end,1,11)){
      cterm = paste0(c(
        substr(library_real_sequence,left_reference_beginning,left_reference_end),
        substr(library_real_sequence,right_reference_beginning,right_reference_end)),collapse = "")
    }
  }else{
    cterm = NA
    cterm_cleavage_pos = NA
  }
  output = data.frame(peptide = peptide_sequence,
                      cterm = cterm,
                      cterm_cleavage_pos = cterm_cleavage_pos)
  return(output)
}
