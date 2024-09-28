#' mspms_tidy
#' Convert a SummarizedExperiment object within a QFeatures object
#' into a tidy tibble.
#'
#' @param processed_qf a QFeature object containing rowData and colData.
#' @param se_name  the name of the SummarizedExperiment you would like to
#' extract
#'
#' @return a tibble containing all the rowData, colData, and assay data for the
#' specified SummarizedExperiment.
#'
#' @export
#' @examples
#' mspms_data <- mspms_tidy(mspms::processed_qf)
mspms_tidy <- function(processed_qf, se_name = "peptides_norm") {
    . <- NULL
    # Extracting assay data
    assay_data <- SummarizedExperiment::assay(processed_qf[[se_name]]) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("peptide") %>%
        tidyr::pivot_longer(2:length(.),
            names_to = "quantCols",
            values_to = se_name
        )
    # Extracting colData
    colData <- SummarizedExperiment::colData(processed_qf) %>%
        SummarizedExperiment::as.data.frame() %>%
        tibble::remove_rownames()
    # Extracting row data
    rowData <- SummarizedExperiment::rowData(processed_qf[[se_name]]) %>%
        SummarizedExperiment::as.data.frame() %>%
        tibble::remove_rownames()
    long <- dplyr::inner_join(colData, assay_data)
    tidy <- dplyr::inner_join(rowData, long, by = "peptide") %>%
        tibble::as_tibble()
    return(tidy)
}
