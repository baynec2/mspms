#' extract_design_matrix_pd
#'
#' extract a design matrix from the proteome discover column headers corresponding to samples.
#'
#' Note this depends on how you have decided to label your samples in the proteome discoverer software.
#' see ../tests/testthat/test-prepare_pd.R for an example.
#'
#' Abundance: F5: Sample, 1, 15 min is an example of a sample name input by the software that this is designed to work with
#'
#' in this example, samples gets assigned to the pd_group, 1 gets assigned to the group, and 15 min gets assigned to the time.
#'
#' @param filepath = this is the path to the proteome discoverer excel output
#'
#' @return a tibble containing the extracted design matrix compatable with the normalyze function
#' @export
#'
#' @examplesIf isTRUE(FALSE)
#'
#' design_matrix = extract_design_matrix_pd("tests/testdata/proteome_discoverer_output.xlsx")

extract_design_matrix_pd = function(filepath){

   # Read the data
  data = readxl::read_excel(filepath)

  # Data col indexes
  col_ind =which(grepl("Abundance.*",names(data)))

  #Find start of data
  start = min(col_ind)
  end = max(col_ind)

  fixed_names = gsub("Abundance: .*: ","",names(data)[start:end])

  # Extract the design matrix
  design_matrix = tibble::tibble(sample = fixed_names) %>%
    tidyr::separate(col = "sample",sep = ",", into = c("pd_group","group","time"), remove = FALSE)

  return(design_matrix)

}
