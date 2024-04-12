#' extract_design_matrix
#'
#'This extracts the design matrix from the column names in the data frame  Note that this assumes a very specific naming convention. 
#'In this naming convention, there are three _. 
#'The first specifies the 
#'The second specifics the time point
#'The third specifies the replicate. 
#'The combination of the first and second defines all of the distinct groups. 
#' @param prepared_data = this is data that has been prepared for normalyzer analysis
#'
#' @return
#' a data frame with the sample (complete names) and group id (samples that are treated as replicates of each other)
#' @export
#'
#' @examples
extract_design_matrix = function(prepared_data){
  
  #Need to extract a design matrix from the sample names, using column names as index
  first_col = which(names(prepared_data) == "Avg. Area") +1
  last_col = which(names(prepared_data) == "Sample Profile (Ratio)") -1
  
  design_matrix = prepared_data %>% 
    # Select the columns that contain the samples
    dplyr::select(first_col:last_col) %>% 
    # put in the long format
    tidyr::pivot_longer(1:length(.)) %>% 
    # spliting sample and group by file name
    tidyr::separate(name,into = c("group1","group2","replicate"),sep = "_") %>% 
    dplyr::mutate(sample = paste0(group1,"_",group2,"_",replicate),group = paste0(group1,"_",group2)) %>% 
    dplyr::select(sample,group) %>% 
    dplyr::distinct() 
  
  return(design_matrix)
  
}