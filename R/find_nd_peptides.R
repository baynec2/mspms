#' find_nd_peptides
#' Find library peptides that were not detected in the data.
#'
#' @param prepared_data a tibble with data prepared by the mspms package.
#' @param peptide_library a tibble containing the peptide library used
#' @param design_matrix a tibble containing the design matrix
#'
#' @return a tibble with the missing peptides for each sample
#' @export
#'
#' @examples
#' nd_peptides <- find_nd_peptides(
#'   mspms::peaks_prepared_data,
#'   mspms::peptide_library,
#'   mspms::design_matrix
#' )
#'
find_nd_peptides <- function(prepared_data, peptide_library, design_matrix) {
  # Which library peptides are never detected?
  data_only <- prepared_data[4:length(prepared_data)]
  `%!in%` <- Negate(`%in%`)

  find_missing <- function(input) {
    peptides_in <- prepared_data$`Protein Accession`[!is.na(input)]
    peptides_out_TF <- peptide_library$library_reference_id %!in% peptides_in
    missing <- peptide_library$library_reference_id[peptides_out_TF]
    out <- tibble::tibble(missing_peptide_id = missing)
    return(out)
  }

  missing_total <- purrr::map_df(data_only, find_missing, .id = "sample") %>%
    dplyr::mutate(type_missing = "completely_missing")

  # Which full length peptides are never detected
  full <- prepared_data %>%
    dplyr::filter(!grepl(".*_.*", Peptide)) 
  
  data_only <- full[4:length(full)]
  
  find_missing_2 <- function(input) {
    peptides_in <- full$`Protein Accession`[!is.na(input)]
    peptides_out_TF <- peptide_library$library_reference_id %!in% peptides_in
    missing <- peptide_library$library_reference_id[peptides_out_TF]
    out <- tibble::tibble(missing_peptide_id = missing)
    return(out)
  }

  missing_full_length <- purrr::map_df(data_only, find_missing_2, .id = "sample") %>%
    dplyr::mutate(type_missing = "undigested_peptide_library")

  out <- dplyr::bind_rows(missing_total, missing_full_length) %>%
    dplyr::inner_join(peptide_library, by = c(
      "missing_peptide_id" =
        "library_reference_id"
    ))

  return(out)
}
