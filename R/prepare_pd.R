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
  # Reading in the label free quantification data
  lfq <- readr::read_delim(peptide_groups_filepath)
  # Check that this file appears to be a proper fragpipe file
  check_file_is_valid_pd(lfq)
  # Reformat the columns to be consistent with mspsms framework
  lfq <- lfq %>%
    dplyr::rename(
      peptide = "Annotated Sequence",
      library_id = "Master Protein Accessions"
    )
  # Converting the peptide notation from proteome discover to standard format.
  lfq <- lfq %>%
    dplyr::mutate(
      peptide = gsub("\\[-\\]", "", .data$peptide),
      peptide = gsub("^\\.", "", .data$peptide),
      peptide = gsub("\\.$", "", .data$peptide),
      peptide = gsub("\\].", "_", .data$peptide),
      peptide = gsub("\\.\\[", "_", .data$peptide),
      peptide = gsub("\\]", "", .data$peptide),
      peptide = gsub("\\[", "", .data$peptide)
    )
  # Converting the column names to reasonable sample ID names
  lfq <- lfq %>%
    dplyr::select(-dplyr::contains("Normalized")) %>%
    dplyr::select(-dplyr::contains("Scaled")) %>% 
    dplyr::select(-dplyr::contains("Count"))
  # Extracting only needed columns
  lfq <- lfq %>%
    dplyr::select("peptide", "library_id", dplyr::contains("Abundance")) %>%
    dplyr::rename_with(~ gsub("Abundance: ", "", .)) %>%
    dplyr::rename_with(~ gsub(": Sample", "", .))
  # Loading colData
  colData <- load_colData(colData_filepath)
  # converting into QF object
  pd_prepared_qf <- prepared_to_qf(
    lfq, colData,
    peptide_library, n_residues
  )
  return(pd_prepared_qf)
}
