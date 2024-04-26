#' extract_design_matrix_pd
#'
#' extract a design matrix from the proteome discover column headers corresponding to samples.
#'
#' @param filepath = this is the path to the proteome discoverer excel output
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#'
#'
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
  design_matrix = tibble(sample = fixed_names) %>%
    tidyr::separate(col = "sample",sep = ",", into = c("pd_group","group","time"), remove = FALSE)

  return(design_matrix)

}
