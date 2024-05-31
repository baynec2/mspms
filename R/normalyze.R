#' normalyze: This function performs normalization of MSP-MS data.
#'
#' This function uses the NormalyzerDE package to perform normalization.
#' It returns the median normalized data as a data frame.
#' For this normalization the Intensity of each variable in a given sample is
#'  divided by the median of intensities of all variables
#' in the sample and then multiplied by the mean of median of sum of
#' intensities of all variables in all samples.
#'
#' @param prepared_data = this is the data that has been prepared via
#' mspms::prepare_for_normalyzer().
#' @param design_matrix = this is the design matrix. This should be a data frame
#'  with a column named "sample" that contains the names of the samples
#'  in the data,
#' a column named "group" that contains the group that the sample belongs to a
#' column named "time" that contains the time the sample was incubated for,
#' and on named "condition" that contains condition cooresponding to each sample
#'  (commonly some kind of protease inhibitor).
#'
#' @param outdir = this is the directory you would like NormalyzerDE::nomralyzer
#' to write output to.Default is current working dir.
#'
#' @return
#' Files from NormalyzerDE are written to the outdir specified.
#'
#' The normalized data is returned as a tibble.
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' normalyzed_data <- normalyze(mspms::peaks_prepared_data, mspms::design_matrix,
#'   outdir = getwd()
#' )
#'
normalyze <- function(prepared_data, design_matrix, outdir = getwd()) {
  # Extracting only the data
  dataOnly <- prepared_data[, design_matrix$sample]

  # Extracting only the annotation columns
  annotOnly <- prepared_data[, !(colnames(prepared_data) %in% design_matrix$sample)]

  # Making a summarized Experiment object
  sumExpObj <- SummarizedExperiment::SummarizedExperiment(
    as.matrix(dataOnly),
    colData = design_matrix,
    rowData = annotOnly
  )


  # Making the jobName be the date folowed by brief tag
  jobName <- paste0(Sys.Date(), "_mspms_normalyze_output")
  ## Performing the normalization, output files to outdir ##
  NormalyzerDE::normalyzer(
    jobName = jobName,
    experimentObj = sumExpObj,
    outdir
  )

  ## Reading in the data frame we just produced that we care about ##

  norm_data <- readr::read_delim(paste0(
    outdir, "/",
    jobName, "/",
    "median-normalized.txt"
  ))

  # The columns after z. correspond to our samples.
  index <- which(names(norm_data) == "Protein Accession") + 1

  # Performing a reverse log2 transformation
  norm_data <- norm_data %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(index:length(.)), ~ 2**.x))

  return(norm_data)
}
