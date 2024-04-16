#' nterm_cleavage
#'
#' Finding the cleavage sequences on the N terminus of a given peptide in reference to the peptide library it was derived from.
#'
#' @param peptide_sequence = this is the peptide sequence in single leter AA code. _ denotes cleavage site..
#' @param library_match_sequence  = this is the sequence of the corresponding peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the Ls transformed to Ms (This is what the legacy code did so it is kept this way in case there was a good reason for it)
#'
#' @return
#' a data frame with the peptide sequence, cleavage sequences 4 AA on the left and right of the n term cleavage, and the position of the n term cleavage in the library sequence.
#' @export
#'
#' @examples
#'
#' peptide_sequence = "F_LAHWVGI"
#' library_real_sequence = "GMYYKRFMAHWVGI"
#' library_match_sequence = "GLYYKRFLAHWVGI"
#'
#'
#' nterm_cleavage(peptide_sequence,library_real_sequence,library_match_sequence)

nterm_cleavage = function(peptide_sequence,library_match_sequence,library_real_sequence){

  # Determining if there is a n term cleavage
  # _ denotes a cleavage, and if it is the second position, it is on the n term!
  if((grepl("_", peptide_sequence)==TRUE) && (regexpr('_',peptide_sequence)[[1]][1]==2)){

    # if there is a n term cleavage, the first letter of the right side is the third letter our sequence according to the logic above.
    pos = 2 + 1

    # taking the sequence from right after the _ to . Will ultimately be a 4 AA long fragment
    temp = substr(peptide_sequence,pos, pos+3)

    # Checking to see what part of the reference sequence this matches.
    right_reference_start = regexpr(temp,library_match_sequence)[[1]][1]
    right_reference_end = right_reference_start + 3

    # the position of the n term clevage is the same as the right reference start - 1
    nterm_cleavage_pos = right_reference_start - 1


    #Now determining the left side of the cleavage event.
    left_reference_end = right_reference_start - 1
    left_reference_start =  left_reference_end - 3

    # Deciding what to do based on what position of the reference the match is to.
    # Note that all of the peptides in the library are 14 AA long
    # This is all about getting the 4 flanking amino acids correct. It should be 4 on each side of the cleavage
    # X denotes that there was nothing there, as happens if it is on the edge of the sequence in the peptide library.

    if(right_reference_start == 2){
      # If the cleavage position is the second AA in the reference sequence, put 3 Xs before hand. Why?
      nterm = paste(
        c("XXX",
          substr(library_real_sequence,left_reference_end,left_reference_end),
          substr(library_real_sequence,right_reference_start,right_reference_end)
          ),collapse = "")
    }else if(right_reference_start == 3){
      nterm = paste(
        c("XX",
          substr(library_real_sequence,left_reference_end-1,left_reference_end),
          substr(library_real_sequence,right_reference_start,right_reference_end)
        ),collapse = "")
    }else if(right_reference_start == 4){
      nterm = paste(
        c("X",
          substr(library_real_sequence,left_reference_end-2,left_reference_end),
          substr(library_real_sequence,right_reference_start,right_reference_end)
        ),collapse = "")
    }else if(dplyr::between(right_reference_start,5,14)){
      nterm = paste(
        c(substr(library_real_sequence,left_reference_start,left_reference_end),
          substr(library_real_sequence,right_reference_start,right_reference_end)
        ),collapse = "")
    }

  }else{
    nterm = NA
    nterm_cleavage_pos = NA
  }
  output = data.frame(peptide = peptide_sequence,
                      nterm = nterm,
                      nterm_cleavage_pos = nterm_cleavage_pos)
  return(output)
  }


