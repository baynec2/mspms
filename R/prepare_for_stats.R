#' prepare_for_stats
#'
#' This function prepares data for downstream statistics by transforming the data to the long format and combining with the design matrix.
#'
#' @param cleavage_added_data = this is the processed data with cleavage information added. Intended to be downstream of the add_cleavages() function.
#' @param design_matrix = this is the design matrix with "sample",group","time","condition" columns.
#'
#' @return a tibble in long format with the design matrix combined.
#' @export
#' @examples
#' prepare_for_stats(mspms::cleavage_added_data,mspms::design_matrix)
prepare_for_stats = function(cleavage_added_data,design_matrix){
  # dealing with no visible binding for global variable ‘.’ NOTE
  . = NULL
  #Figuring out where the samples start
  start_of_samples = which(names(cleavage_added_data) == "cterm_cleavage_pos")+1

  # make into long format, combine with design matrix.
  long = cleavage_added_data %>%
    tidyr::pivot_longer(start_of_samples:length(.),names_to = "sample") %>%
    dplyr::inner_join(design_matrix,by = "sample") %>%
    tibble::as_tibble()

  return(long)

}
