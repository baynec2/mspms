#' prepare_sage
#' Prepare a label free quantification file exported from Sage
#' for subsequent mspms analysis.
#' @param sage_lfq_filepath filepath to lfq.tsv file output from
# ` Sage.
#' @param colData_filepath file path to .csv file containing colData.
#' Must have columns named "quantCols","group","condition",and "time".
#' @param peptide_library peptide library used with experiment. Contains
#' columns "library_id", "library_match_sequence", and "library_real_sequence".
#' @param n_residues the number of amino acid residues before and after the
#' cleavage site to generate a cleavage seq for.
#' @return a QFeatures object containing a summarizedExperiment named "peptides"
#' @export
#' @examples
#' sage_lfq_filepath <- system.file(
#'   "extdata/sage_lfq.tsv",
#'   package = "mspms"
#' )
#' colData_filepath <- system.file("extdata/colData.csv", package = "mspms")
#'
#' prepare_sage(sage_lfq_filepath, colData_filepath)
prepare_sage <- function(sage_lfq_filepath,
                         colData_filepath,
                         peptide_library = mspms::peptide_library,
                         n_residues = 4) {
  prepared_data <- prepare_file(
    filepath = sage_lfq_filepath,
    colData_filepath = colData_filepath,
    peptide_library = peptide_library,
    n_residues = n_residues,
    read_fun = readr::read_tsv,
    validate_fun = check_file_is_valid_sage,
    transform_fun = transform_sage
  )

  return(prepared_data)
}
