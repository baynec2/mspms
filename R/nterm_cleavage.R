#' nterm_cleavage
#'
#' Finding the cleavage sequences on the N terminus of a given peptide in reference to the peptide library it was derived from.
#'
#' @param peptide_sequence = this is the peptide sequence in single leter AA code. _ denotes cleavage site..
#' @param library_match_sequence  = this is the sequence of the corresponding peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the Ls transformed to Ms (This is what the legacy code did so it is kept this way in case there was a good reason for it)
#' @param n_residues = this is the number of residues to the left and right of the cleavage site that you want to extract.

#' @return a tibble with the peptide sequence, cleavage sequences n specified number of AA on the left and right of the n term cleavage, and the position of the n term cleavage in the library sequence.
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

nterm_cleavage =function(peptide_sequence,library_match_sequence,library_real_sequence,n_residues = 4){

  # Determining if there is a n term cleavage
  # _ denotes a cleavage, and if it is the second position, it is on the n term!
  if((grepl("_", peptide_sequence)==TRUE) && (regexpr('_',peptide_sequence)[[1]][1]==2)){

    # if there is a n term cleavage, the first letter of the right side is the third letter our sequence according to the logic above.
    pos = 2 + 1

    # taking the sequence from right after the _ to . Will ultimately be a 4 AA long fragment
    temp = substr(peptide_sequence,pos, pos+(n_residues -1))

    # Checking to see what part of the reference sequence this matches.
    right_reference_start = regexpr(temp,library_match_sequence)[[1]][1]
    right_reference_end = right_reference_start + (n_residues - 1)

    # the position of the n term clevage is the same as the right reference start - 1
    nterm_cleavage_pos = right_reference_start - 1

    #Now determining the left side of the cleavage event.
    left_reference_end = right_reference_start - 1
    left_reference_start =  left_reference_end - (n_residues - 1)

    # Extracting the sequences from the reference sequence
    right_sequence = substr(library_real_sequence,right_reference_start,right_reference_end)
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    right_sequence = paste0(right_sequence,paste0(rep("X",n_residues - nchar(right_sequence),collapse = "")))

    left_sequence = substr(library_real_sequence,left_reference_start,left_reference_end)
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    left_sequence = paste0(paste0(rep("X",n_residues - nchar(left_sequence)),collapse = ""),left_sequence)


    nterm = paste(c(left_sequence,right_sequence),collapse = "")

  }else{
    nterm = NA
    nterm_cleavage_pos = NA
  }
  output = tibble::tibble(peptide = peptide_sequence,
                      nterm = nterm,
                      nterm_cleavage_pos = nterm_cleavage_pos)
  return(output)
}


