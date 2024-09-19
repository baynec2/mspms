#' Find library peptides that were not detected in the data.
#'
#' @param prepared_data a tibble with data prepared by the mspms package.
#' @param peptide_library a tibble containing the peptide library used
#' @param design_matrix a tibble containing the design matrix
#'
#' @return a tibble with the missing peptides for each sample
#' @keywords internal
#'
#' @examples
#' nd_peptides <- find_nd_peptides(
#'   mspms::peaks_prepared_data,
#'   mspms::peptide_library,
#'   mspms::design_matrix
#' )
#'
calc_per_samples_library_nd <- function(processed_qf,
                                        peptide_library_ids = mspms::peptide_library$library_id) {
  # tidying the data so we can work with it
  mspms_data = mspms_tidy(processed_qf,"peptides")
  #calculating n of samples
  n_samples = length(unique(mspms_data$quantCols))
  # converting to wide format, easier to reason with this way
  wide = mspms_data %>% 
    dplyr::select("peptide","library_id","peptides",
                  "peptide_type","quantCols") 
  # considering full length peptides only   
  full_length = wide %>% 
    dplyr::filter(peptide_type == "full_length") %>% 
    tidyr::pivot_wider(names_from = "quantCols", 
                       values_from = "peptides")
  # What peptides are not detected at all
  nd_full = peptide_library_ids[peptide_library_ids %!in% unique(
    full_length$library_id)]
  # Building a tibble of the peptides not detected
  nd_full = tibble::tibble(library_id = nd_full,
                           n_samples = n_samples,
                           n_missing = n_samples)
# Figuring the number of full length library ids missing per sample
  n_missing_full_length = tibble::tibble(library_id = full_length$library_id,
                                         n_samples = n_samples,
                                         n_missing = rowSums(is.na(
                                           full_length)))
  # Combining the completely missing with partially missing data
  full = dplyr::bind_rows(n_missing_full_length,nd_full) %>% 
    dplyr::mutate(peptide_type ="full_length",.after = "library_id")
  # Now considering cleavage products, conting the number of non nas, per sample
  cleavage_product = wide %>% 
    dplyr::filter(peptide_type == "cleavage_product") %>% 
    dplyr::select(-"peptide") %>% 
    tidyr::pivot_wider(names_from = quantCols, 
                       values_from = peptides,
                       values_fn = ~sum(!is.na(.)))
  # If there are 0 non NAs, there must only be NA values for that sample
  cleavage_product[cleavage_product == 0] <- NA
  # Creating a tibble with all the data
  cp_row_sums = tibble::tibble(library_id = cleavage_product$library_id,
                               n_samples = n_samples,
                               n_missing = rowSums(is.na(cleavage_product)))
  
  nd_cleavage = peptide_library_ids[peptide_library_ids %!in% unique(
    cleavage_product$library_id)]
  
  nd_cleavage = tibble::tibble(library_id = nd_cleavage,
                               n_samples = n_samples,
                               n_missing = n_samples)
  cp_final = dplyr::bind_rows(cp_row_sums,
                              nd_cleavage) %>% 
    dplyr::mutate(peptide_type = "cleavage_product",.after = "library_id")
  #Combining all data 
  out = dplyr::bind_rows(full,cp_final) %>% 
    dplyr::mutate(per_samples_undetected = n_missing/ n_samples * 100)
  return(out)
}
#' prepare_qc_check
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
#' 
#' @keywords internal
#'
#' @examples
#' qc_check_data <- qc_check(
#'   mspms::peaks_prepared_data,
#'   mspms::peptide_library,
#'   mspms::design_matrix
#' )
prepare_qc_check_data <- function(processed_qf,
                                  peptide_library = mspms::peptide_library) {
  
  # Filtering to only include detected peptides
  long <- processed_qf %>%
    mspms_tidy(se_name = "peptides") %>% 
    dplyr::filter(!is.na(.data$peptides))
  
  # Counting # of peptides in peptide library
  library_num <- length(peptide_library$library_id)
  
  # checking to see what percentage of the library is detected in each sample.
  check <- long %>%
    dplyr::group_by(.data$quantCols,
                    .data$peptide_type,
                    .data$condition,
                    .data$time,
                    .data$group) %>%
    dplyr::summarise(n_detected = n_distinct(.data$library_id),
                     n_total = library_num,
                     per_library_id_detected = round(n_detected / n_total
                                                     * 100,2),
                     per_library_id_undetected = 100 - 
                       per_library_id_detected)
  
  return(check)
}
