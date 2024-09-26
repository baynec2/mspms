#' all_possible_8mers_from_228_library
#' All possible 8mers from the standard (as of 26April2024)
#' 228 MSP-MS peptide library
#' (This is equivalent to the result of
#' mspms::calculate_all_cleavages(mspms::peptide_library$real_cleavage_seq,n=4))
#' vector of the 14 AA peptides used in the library.
#' @format ## `all_possible_8mers_from_228_library`
#' A vector with 2964 entries
#' @source <standard peptide library used with MSP-MS method
#' in the O’Donoghue lab as of 26April2024>
"all_possible_8mers_from_228_library"

#' peptide_library
#'
#' This is the 228 peptide library used by the O’Donoghue lab as of 26April2024.
#' @format ## `peptide_library`
#' A  data frame with 228 rows and 3 columns:
#' \describe{
#'   \item{library_reference_id}{reference id of the detected peptide as
#'   put in upstream software}
#'   \item{library_match_sequence}{the sequence match to the peptide library,
#'    methionine is replaced with norleucine,which should function the same as
#'     methionine for proteases but has the same mass as L}
#'   \item{library_real_sequence}{Ls corresponding to norleucine are replaced
#'    back with n (for norleucine )}
#'   ...
#' }
#' @source <O’Donoghue lab as of 26April2024 >
"peptide_library"
#' peaks_prepared_data
#' A QFeatures object prepared from PEAKS data of cathepsin data/.
#' @format ## `peaks_prepared_data`
#' An instance of class QFeatures containing 1 assays:
#' [1] peptides: SummarizedExperiment with 2071 rows and 42 columns
#' \describe{
#'   \item{peptides}{Peptide Sequence Detected}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in
#'  "tests/testdata/protein-peptides-id.csv" and
#'  "tests/testdata/protein-peptides-lfq.csv">
"peaks_prepared_data"

#' processed_qf
#' A QFeatures object prepared from PEAKS data of Cathepsin data that has been
#' processed (imputation/normalization)
#' @format ## `peaks_prepared_data`
#' An instance of class QFeatures containing 5 assays:
#' [1] peptides: SummarizedExperiment with 2071 rows and 42 columns
#' [2] peptides_log: SummarizedExperiment with 2071 rows and 42 columns
#' [3] peptides_log_norm: SummarizedExperiment with 2071 rows and 42 columns
#' [4] peptides_log_impute_norm: SummarizedExperiment with 2071 rows and 42
#' columns
#' [5] peptides_norm: SummarizedExperiment with 2071 rows and 42 columns
#' \describe{
#'   \item{peptides}{Peptide Sequence Detected}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in
#'  "tests/testdata/protein-peptides-id.csv" and
#'  "tests/testdata/protein-peptides-lfq.csv">
"processed_qf"


#' log2fc_t_test_data
#' A tibble containing the results of t-tests and log2fc compared to time 0
#' 14,497 × 19
#' @format ## `peaks_prepared_data`
#' A tibble: 14,497 × 19
#' @source <mspms processed data originally from PEAKS files found in
#'  "tests/testdata/protein-peptides-id.csv" and
#'  "tests/testdata/protein-peptides-lfq.csv">
"log2fc_t_test_data"

#' colData
#' A tibble containing the colData associated with an experiment to proc
#' @format ## `colData`
#' A tibble: 42 × 4
#' @source colData corresponding to cathepsin A-D MSP-MS experiment
"colData"

#' mspms_tidy_data
#' A tibble containing tidy data derived from QFeatures object
#' @format ## `mspms_tidy_data`
#' A tibble:
#' @source processed_qf
"mspms_tidy_data"
