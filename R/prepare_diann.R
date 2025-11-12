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
#' colData_filepath <- system.file("extdata/diann_colData.csv", package = "mspms")
#' prepare_diann(precursor_filepath, colData_filepath)
prepare_diann <- function(precursor_filepath,
                          colData_filepath,
                          peptide_library = mspms::peptide_library,
                          n_residues = 4) {
  prepared_data <- prepare_file(
    filepath = precursor_filepath,
    colData_filepath = colData_filepath,
    peptide_library = peptide_library,
    n_residues = n_residues,
    read_fun = readr::read_tsv,
    validate_fun = check_file_is_valid_diann,
    transform_fun = transform_diann
  )

  return(prepared_data)
}
