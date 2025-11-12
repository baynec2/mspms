#' Prepare PEAKS label-free quantification data for MSP-MS analysis
#'
#' This function reads, validates, transforms, and converts a PEAKS
#' LFQ file into a `QFeatures` object compatible with the `mspms` workflow.
#'
#' @param lfq_filepath Path to the PEAKS `.csv` file containing peptide-level LFQ data.
#' @param colData_filepath Path to a `.csv` file containing sample metadata (`colData`).
#'   Must include the columns `"quantCols"`, `"group"`, `"condition"`, and `"time"`.
#' @param quality_threshold Minimum quality score required for a peptide to be retained.
#'   Peptides below this threshold are filtered out (default `0.3`).
#' @param peptide_library A peptide library used in the experiment, typically
#'   `mspms::peptide_library`. Must include `"library_id"`, `"library_match_sequence"`,
#'   and `"library_real_sequence"`.
#' @param n_residues Number of amino acid residues to include on each side of the
#'   cleavage site when generating cleavage sequences (default `4`).
#'
#' @return A `QFeatures` object containing a `SummarizedExperiment` named `"peptides"`.
#' @export
#'
#' @examples
#' lfq_filepath <- system.file(
#'   "extdata/peaks_protein-peptides-lfq.csv",
#'   package = "mspms"
#' )
#' colData_filepath <- system.file(
#'   "extdata/colData.csv",
#'   package = "mspms"
#' )
#' peaks_qf <- mspms::prepare_peaks(lfq_filepath, colData_filepath)
prepare_peaks <- function(lfq_filepath,
                          colData_filepath,
                          quality_threshold = 0.3,
                          peptide_library = mspms::peptide_library,
                          n_residues = 4) {
  # Prepare PEAKS data using the generic prepare_file workflow
  prepared_data <- prepare_file(
    filepath = lfq_filepath,
    colData_filepath = colData_filepath,
    peptide_library = peptide_library,
    n_residues = n_residues,
    read_fun = readr::read_csv,
    validate_fun = check_file_is_valid_peaks,
    transform_fun = function(df, peptide_library) {
      transform_peaks(df, peptide_library, quality_threshold)
    }
  )

  return(prepared_data)
}
