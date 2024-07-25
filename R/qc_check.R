#' qc_check
#' Run simple quality control checks on the data. This checks to see how many
#' full length peptides belonging to the library were identified in the data,
#' as well as how many peptides were detected among all cleavage products
#'
#' @param prepared_data data prepared for analysis by the mspms package.
#' @param peptide_library the peptide library used in the experiment
#' @param design_matrix the design matrix containing experimental metadata
#' @return a list with two slots, 1st slot contains a tibble with QC metrics
#' second list contains and interpretation. 2nd slot contains a tibble with the
#' missing peptides in each sample.
#' @export
#'
#' @examples
#' qc_check_data <- qc_check(
#'   mspms::peaks_prepared_data,
#'   mspms::peptide_library,
#'   mspms::design_matrix
#' )
qc_check <- function(prepared_data, peptide_library, design_matrix) {
  # Preparing the data for QC analysis
  long <- prepared_data %>%
    # column 4 is the start of the data
    tidyr::pivot_longer(4:length(.), names_to = "sample") %>%
    dplyr::filter(!is.na(value))
  # Counting # of peptides in library
  library_num <- length(peptide_library$library_reference_id)
  # checking to see how many full length peptides were found in each sample
  check_1 <- long %>%
    dplyr::filter(!grepl(".*_.*", Peptide)) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(n_undigested_library_detected = dplyr::n_distinct(Peptide))
  # Checking to see how many components of the library are detected total
  check_2 <- long %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_library_detected =
        dplyr::n_distinct(`Protein Accession`)
    )
  # Combining the checks
  out <- dplyr::right_join(check_1, check_2, by = "sample") %>%
    dplyr::mutate(
      perc_undigested_library_detected =
        round((n_undigested_library_detected / library_num * 100), 1),
      perc_library_detected = round(
        n_library_detected / library_num * 100, 1
      )
    ) %>%
    dplyr::inner_join(design_matrix, by = "sample")

  return(out)
}
