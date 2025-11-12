#' prepare_pd
#' Prepare a label free quantification file exported from Proteome Discoverer
#' for subsequent mspms analysis.
#' @param peptide_groups_filepath filepath to PeptideGroups.txt file exported
#' from proteome discoverer.
#' @param colData_filepath file path to .csv file containing colData.
#' Must have columns named "quantCols","group","condition",and "time".
#' @param peptide_library peptide library used with experiment. Contains
#' columns "library_id", "library_match_sequence", and "library_real_sequence".
#' @param n_residues the number of amino acid residues before and after the
#' cleavage site to generate a cleavage seq for.
#' @return a QFeatures object containing a summarizedExperiment named "peptides"
#' @export
#' @examples
#' peptide_groups_filepath <- system.file(
#'   "extdata/proteome_discoverer_PeptideGroups.txt",
#'   package = "mspms"
#' )
#' colData_filepath <- system.file("extdata/colData.csv", package = "mspms")
prepare_pd <- function(peptide_groups_filepath,
                       colData_filepath,
                       peptide_library = mspms::peptide_library,
                       n_residues = 4) {
  prepared_data <- prepare_file(
    filepath = peptide_groups_filepath,
    colData_filepath = colData_filepath,
    peptide_library = peptide_library,
    n_residues = n_residues,
    read_fun = readr::read_delim,
    validate_fun = check_file_is_valid_pd,
    transform_fun = transform_pd
  )

  return(prepared_data)
}
