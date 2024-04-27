#' cterm_cleavage
#'
#' Finding the cleavage sequences on the C terminus of a given peptide in reference to the peptide library it was derived from
#'
#' @param peptide_sequence = this is the peptide sequence in single leter AA code. _ denotes cleavage site.
#' @param library_match_sequence  = this is the sequence of the corresponding peptide in the library to the peptide_sequence.
#' @param library_real_sequence = this is the library match sequence with the Ls transformed to Ms (This is what the legacy code did so it is kept this way in case there was a good reason for it)
#' @param n_residues = this is the number of residues to the left and right of the cleavage site that you want to extract.

#' @return a tibble with the peptide sequence, cleavage sequences n number of AA on the left and right of the c term cleavage, and the position of the c term cleavage in the library sequence.
#' @export
#'
#' @examples
#' peptide_sequence = "YWLSTHLAGKR_R"
#' library_match_sequence = "YWLSTHLAGKRRDW"
#' library_real_sequence = "YWMSTHLAGKRRDW"
#' cterm_cleavage(peptide_sequence,library_match_sequence,library_real_sequence)
#'

cterm_cleavage = function(peptide_sequence,library_match_sequence,library_real_sequence,n_residues = 4){

  # First, we determine if there is a cterm cleavage
  # We do this by looking to see if there is a _ (denotes a cleavage)
  # If there is, we check to see if it is in the second to last position of the peptide when there is only one _ (only cterm cleavage)
  # Or there might be 2 _s in the peptide sequence (meaning there is a cterm and nterm cleavage)

  n_of_match = length(gregexpr('_',peptide_sequence)[[1]])
  if((grepl("_", peptide_sequence)==TRUE) && (gregexpr('_',peptide_sequence)[[1]][n_of_match]==nchar(peptide_sequence)-1)){

    # if there is a c term cleavage, it is the last position - 1 of the peptide sequence
    pos = nchar(peptide_sequence) - 1

    # now we define the sequence on the left side of the cleavage
    temp = substr(peptide_sequence,pos-n_residues, pos-1)

    # Checking to see what part of the reference sequence this matches.
    left_reference_beginning = regexpr(temp,library_match_sequence)[[1]][1]
    # it is four AA long, so the end is the beginning + 3
    left_reference_end = left_reference_beginning + (n_residues - 1)

    # now we can determine the right side of the cleavage sequence by taking the left reference end and adding 1
    right_reference_beginning = left_reference_end + 1

    # the cterm cleavage position is the right reference beginning minus 1
    cterm_cleavage_pos = right_reference_beginning -1

    # it is four AA long, so the end is the beginning + 3
    right_reference_end = right_reference_beginning + (n_residues - 1)

    # Extracting the sequences from the reference sequence
    right_sequence = substr(library_real_sequence,right_reference_beginning,right_reference_end)
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    right_sequence = paste0(right_sequence,paste0(rep("X",n_residues - nchar(right_sequence)),collapse = ""))

    # Doing the same on the left side of the clevage event

    left_sequence = substr(library_real_sequence,left_reference_beginning,left_reference_end)
    # now we need to make sure that we add X to represent cases where there was no more sequences in the library peptide
    left_sequence = paste0(paste0(rep("X",n_residues - nchar(left_sequence)),collapse = ""),left_sequence)


    cterm = paste(c(left_sequence,right_sequence),collapse = "")
  }else{
    cterm = NA
    cterm_cleavage_pos = NA
  }
  output = tibble::tibble(peptide = peptide_sequence,
                      cterm = cterm,
                      cterm_cleavage_pos = cterm_cleavage_pos)
  return(output)
}
