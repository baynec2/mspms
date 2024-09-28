#' prepare_peaks
#' Prepare a label free quantification file exported from PEAKS
#' for subsequent mspms analysis.
#' @param lfq_filepath  this is the file path to a .csv file exported from
#' PEAKS
#' @param colData_filepath file path to .csv file containing colData.
#' Must have columns named "quantCols","group","condition",and "time".
#' @param quality_threshold only consider peptides with quality scores > than
#' this threshold.
#' @param peptide_library peptide library used in the experiment.
#' @param n_residues the number of amino acid residues before and after the
#' cleavage site to generate a cleavage seq for.
#' @return a QFeatures object containing a summarizedExperiment named "peptides"
#' @export
#' @examples
#' lfq_filepath <- system.file("extdata/peaks_protein-peptides-lfq.csv", package = "mspms")
#' colData_filepath <- system.file("extdata/colData.csv", package = "mspms")
#' # Prepare the data
#' peaks_prepared_data <- mspms::prepare_peaks(lfq_filepath, colData_filepath)
prepare_peaks <- function(lfq_filepath,
                          colData_filepath,
                          quality_threshold = 0.3,
                          peptide_library = mspms::peptide_library,
                          n_residues = 4) {
    # Reading in the label free quantification data
    lfq <- readr::read_csv(lfq_filepath, na = c("-", 0, "", "NA"))
    # Check that this file appears to be a proper peaks file
    check_file_is_valid_peaks(lfq)
    # Filtering out low quality peptides
    quality <- lfq %>%
        dplyr::filter(.data$Quality > quality_threshold)
    # Letting the user know how many peptides were filtered out with threshold
    n_peptides_rm <- sum(lfq$Quality < quality_threshold, na.rm = TRUE)
    message <- " peptides were removed because they had a quality score < "
    percentage <- round(n_peptides_rm / length(lfq$Peptide) * 100, 0)
    # Printing how many peptides are filtered out
    print(paste0(n_peptides_rm, message, quality_threshold, " (", percentage, "%)"))
    # Selecting columns containing data
    start <- which(names(lfq) == "Avg. Area") + 1
    end <- which(names(lfq) == "Sample Profile (Ratio)") - 1
    # Only reporting relevant columns
    selected <- lfq %>%
        dplyr::mutate(Peptide = gsub("\\.", "_", .data$Peptide)) %>%
        dplyr::select(
            "peptide" = "Peptide", "library_id" = "Protein Accession",
            dplyr::any_of(start:end)
        )
    # reading in colData
    colData <- load_colData(colData_filepath)
    # converting into QF object
    peaks_prepared_qf <- prepared_to_qf(
        selected, colData,
        peptide_library, n_residues
    )
    return(peaks_prepared_qf)
}
