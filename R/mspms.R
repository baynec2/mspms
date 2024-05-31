#' mspms
#' runs the standard mspms data analysis pipeline in it's entirety.
#'
#' @param prepared_data = this is data that has been processed by one of the prepare data functions (which one depends on the upstream data processing software used.)
#' @param design_matrix = this is the design matrix. This should be a data frame with a column named "sample" that contains the names of the samples in the data,
#' a column named "group" that contains the group that the sample belongs to a column named "time" that contains the time the sample was incubated for,
#' and on named "condition" that contains condition cooresponding to each sample (commonly some kind of protease inhibitor).
#' @param peptide_library = this is the peptide library that has been used in the experiment.
#' @param n_residues = this is the number of cleavages to the left and right of the cleavage site you want to include in the cleavage motif.
#' @param outdir = this is the directory to write the results of normalyze to.
#' @return a tibble with processed mspms data.
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' out = mspms::mspms(mspms::peaks_prepared_data,mspms::design_matrix)
#'
mspms = function(prepared_data,
                 design_matrix,
                 peptide_library = mspms::peptide_library,
                 n_residues = 4,
                 outdir = getwd()){
  #Normalyzing data
  normalyzed_data = mspms::normalyze(prepared_data,design_matrix,outdir)
  #handle outliers, impute
  outliers = mspms::handle_outliers(normalyzed_data,design_matrix)
  imputed = mspms::impute(outliers)
  # Join with the library
  joined_with_library = mspms::join_with_library(imputed,)
  # Add cleavage information
  cleavage_added_data = mspms::add_cleavages(joined_with_library,n_residues = 4)

  #Prepare for stats
  prepared_for_stats = mspms::prepare_for_stats(cleavage_added_data,design_matrix)

  # Polish for downstream analysis (gets rid of cases where there are both n and c cleavage)
  polished_data = mspms::polish(prepared_for_stats)

  return(polished_data)
}
