#' generate_report
#'
#' wrapper function to generate an automatic .html report of a basic mspms
#' analysis.
#'
#' @param prepared_data a QFeatures object containing a SummarizedExperiment
#' named "peptides".
#' @param peptide_library peptide library used with experiment. Contains
#' columns "library_id", "library_match_sequence", and "library_real_sequence".
#' @param n_residues the number of amino acid residues before and after the
#' cleavage site to generate a cleavage seq for.
#' @param outdir the output directory you would like to render the
#' report to.
#' @param output_file the file name to export.
#'
#' @return a knited .html report of the mspms analysis.
#' @export
#' @examplesIf isTRUE(FALSE)
#' generate_report(mspms::peaks_prepared_data)
generate_report <- function(prepared_data,
                            peptide_library = mspms::peptide_library,
                            n_residues = 4,
                            outdir = getwd(),
                            output_file = paste0(
                                Sys.Date(),
                                "_mspms_report.html"
                            )) {
    rmarkdown::render(
        system.file("rmarkdown/templates/mspms_report/skeleton/skeleton.RMD",
            package = "mspms"
        ),
        params = list(
            prepared_data = prepared_data,
            peptide_library = peptide_library,
            n_residues = n_residues
        ),
        clean = TRUE,
        output_dir = outdir,
        output_file = output_file
    )
}
