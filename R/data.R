#' all possible 8mers from the standard 228 peptide library that can
#'
#' equivalent to the result of callling mspms::calculate_all_cleavages() on a
#'  vector of the 14 AA peptides used in the library.
#' @format ## `all_possible_8mers_from_228_library`
#' A vector with 2964 entries
#' @source <Peptide library used in the O’Donoghue lab as of 26April2024>
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
#' this is an example dataset from peaks that has been processed by mspms.
#' @format ## `peaks_prepared_data`
#' A tibble with 820 rows and 27 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in
#'  "tests/testdata/protein-peptides-id.csv" and
#'  "tests/testdata/protein-peptides-lfq.csv">
"peaks_prepared_data"



