#' data prepared for statistics by the mspms package
#'
#' data prepared for statistics by the mspms package as described in the readme
#' @format ## `prepared_for_stats`
#' A tibble with 19,680 rows and 13 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   \item{library_reference_id}{reference id of the detected peptide to that of the library reported by software}
#'   \item{library_match_sequence}{the sequence match to the peptide library}
#'   \item{library_real_sequence}{the real sequence of the peptide in the library}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"prepared_for_stats"

#' all possible 8mers from the 228 peptide library used by the O’Donoghue lab as of 26April2024.
#'
#' equivalent to the result of callling mspms::calculate_all_cleavages() on a vector of the 14 AA peptides used in the library.
#' @format ## `all_possible_8mers_from_228_library`
#' A vector with 2964 entries
#' @source <Peptide library used in the O’Donoghue lab as of 26April2024>
"all_possible_8mers_from_228_library"


#' data with cleavage information added by the mspms package
#'
#' data prepared by the mspms package with cleavage information added as described in the readme
#' @format ## `cleavage_added_data`
#' A tibble with 820 rows and 32 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   \item{library_reference_id}{reference id of the detected peptide to that of the library reported by software}
#'   \item{library_match_sequence}{the sequence match to the peptide library}
#'   \item{library_real_sequence}{the real sequence of the peptide in the library}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"cleavage_added_data"


#' design matrix
#'
#' example design matrix designed to be used with the peaks processed data in the readme
#' @format ## `design_matrix`
#' A tibble with 24 rows and 4 columns:
#' \describe{
#'   \item{sample}{sample name as reported by upstream software }
#'   \item{group}{information about what group the sample belongs to, should be the same for all replicates}
#'   \item{condition}{the condition the sample was exposed to (commonly a protease inhibitor treatment for example)}
#'   \item{time}{the amount of time the sample was incubated for}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"design_matrix"



#' imputed
#'
#' data imputed by the mspms package as described in the readme
#' @format ## `imputed`
#' A tibble with 820 rows and 28 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   \item{RT}{Retention Time}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"imputed"


#' normalyzed_data
#'
#' data normalyzed by the mspms package as described in the readme
#' @format ## `normalyzed_data`
#' A tibble with 820 rows and 27 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   \item{RT}{Retention Time}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"normalyzed_data"


#' outliers
#'
#' data with outliers removed by the mspms package as described in the readme
#' @format ## `outliers`
#' A tibble with 19,680 rows and 5 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   \item{RT}{Retention Time}
#'   \item{Protein Accession}{What each peptide was detected as belonging to by the upstream data processing software}
#'   \item{sample_id}{the sample identifier}
#'   \item{value}{normalyzed values}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"outliers"


#' peptide_library
#'
#' This is the 228 peptide library used by the O’Donoghue lab as of 26April2024.
#' @format ## `peptide_library`
#' A  data frame with 228 rows and 3 columns:
#' \describe{
#'   \item{library_reference_id}{reference id of the detected peptide as put in upstream software}
#'   \item{library_match_sequence}{the sequence match to the peptide library, methionine is replaced with norleucine,which should function the same as methionine for proteases but has the same mass as L}
#'   \item{library_real_sequence}{Ls corresponding to norleucine are replaced back with methionine (even though it is norleucine in reality)}
#'   ...
#' }
#' @source <O’Donoghue lab as of 26April2024 >
"peptide_library"


#' peaks_prepared_data
#' this is an example dataset from peaksthat has been processed by mspms and is ready for normalyzation.
#' @format ## `peaks_prepared_data`
#' A tibble with 820 rows and 27 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"peaks_prepared_data"


#' joined_with_library
#' this is an example dataset processed by mspms that has been joined with the peptide library
#' @format ## `joined_with_library`
#' A tibble with 820 rows and 30 columns:
#' \describe{
#'   \item{Peptide}{Peptide Sequence Detected}
#'   ...
#' }
#' @source <mspms processed data originally from PEAKS files found in "tests/testdata/protein-peptides-id.csv" and "tests/testdata/protein-peptides-lfq.csv">
"joined_with_library"




