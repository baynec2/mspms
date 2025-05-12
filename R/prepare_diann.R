#' prepare_diann
#' 
#' prepare data from the pr_matrix.tsv diann output. This can be either from
#' DIA-NN or from Fragpipe (as it uses DIA-NN for quantification internally for
#' MSFragger-DIA workflows)
#'
#' @param precursor_filepath filepath to report.pr_matrix.tsv file exported
#' from DIA-NN.
#' @param colData_filepath file path to .csv file containing colData.
#' Must have columns named "quantCols","group","condition",and "time".
#' @param peptide_library peptide library used with experiment. Contains
#' columns "library_id", "library_match_sequence", and "library_real_sequence".
#' @param n_residues the number of amino acid residues before and after the
#' cleavage site to generate a cleavage seq for.
#'
#' @returns a QFeatures object. 
#' @export
#'
#' @examples
#' precursor_filepath <- system.file(
#'   "extdata/diann_report.pr_matrix.tsv",
#'   package = "mspms"
#' )
#' colData_filepath <- system.file("extdata/diann_colData.csv",package = "mspms")
#' prepare_diann(precursor_filepath,colData_filepath)
prepare_diann <- function(precursor_filepath,
                          colData_filepath,
                          peptide_library = mspms::peptide_library,
                          n_residues = 4){
  
  # Reading in the quantification data
  lfq <- readr::read_tsv(precursor_filepath)
  
  # Check that this file appears to be a proper diann file.
  check_file_is_valid_diann(lfq)
  
  # start and and end indicies of quantification columns
  start = which(names(lfq) == "Precursor.Id") + 1
  end = length(names(lfq)) -2
  nsamples = end-start
  
  # Convert to peptides by summarizing precursor intensities (sum) per peptide
  peptides <- lfq %>%
    dplyr::group_by(
      library_id = `Protein.Group`,
      peptide = `Stripped.Sequence`
    ) %>%
    dplyr::summarise(
      dplyr::across(start:end, sum, na.rm = TRUE),
      .groups = "drop"
    ) 
  
  peptides <- dplyr::inner_join(peptide_library, peptides, by = "library_id") %>%
    dplyr::mutate(
      cleavage_pos = stringr::str_locate(library_match_sequence, peptide),
      prev_aa_i = cleavage_pos[, "start"] - 1,
      next_aa_i = cleavage_pos[, "end"] + 1,
      prev_aa = dplyr::if_else(prev_aa_i > 0, stringr::str_sub(library_match_sequence, prev_aa_i, prev_aa_i), ""),
      next_aa = dplyr::if_else(next_aa_i <= stringr::str_length(library_match_sequence), stringr::str_sub(library_match_sequence, next_aa_i, next_aa_i), "")
    ) %>% 
    dplyr::mutate(peptide = dplyr::case_when(prev_aa != "" ~ paste0(prev_aa,"_",peptide),
                                             next_aa != "" ~ paste0(peptide,"_",next_aa),
                                             TRUE ~ peptide))
  peptides = peptides %>% 
    dplyr::select(1:(length(names(peptides)) - 5)) %>% 
    dplyr::select(peptide,library_id,dplyr::everything()) %>% 
    dplyr::select(-library_match_sequence,-library_real_sequence)
  
  # Loading colData
  colData <- load_colData(colData_filepath)
  # converting into QF object
  diann_prepared_qf <- prepared_to_qf(
    peptides, colData,
    peptide_library, n_residues
  )
  return(diann_prepared_qf)
}
