#' generate_report
#'
#' wrapper function to generate an automatic .html report of a basic mspms analysis.
#'
#' @param prepared_data = this is the prepared data.
#' @param design_matrix = this is the design matrix.
#' @param peptide_library = this is the peptide library used in the experiment.
#' @param n_residues = this is the number of residues that you want to show cleavage sequences for left and right of the cut site.
#' @param outdir = this is the output directory you would like to render the report to.
#'
#' @return a knited .html report of the mspms analysis.
#' @export
#'
#' @examplesIf is.TRUE(FALSE)
#'
#' generate_report(prepared_data = mspms::peaks_prepared_data,
#' design_matrix = mspms::design_matrix,
#' outdir = "../Desktop/test_report")

generate_report = function(prepared_data,
                           design_matrix,
                           peptide_library = mspms::peptide_library,
                           n_residues = 4,
                           outdir = getwd()){

  rmarkdown::render("inst/rmarkdown/templates/mspms_report/skeleton/skeleton.Rmd",
                    params = list(prepared_data = prepared_data,
                                  design_matrix = design_matrix,
                                  peptide_library = peptide_library,
                                  n_residues = n_residues),
                    clean = TRUE,
                    output_dir = outdir,
                    output_file = paste0(Sys.Date(),"_mspms_report.html"))

}



