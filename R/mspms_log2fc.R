#' mspms_log2fc
#'
#' Calculates the log2fc for each timepoint within each condition relative to time 0.
#'
#' @param prepared_for_stats = this is the data that has been prepared for stats using the prepare for stats function.
#'
#' @return a data frame with the log2fc for each timepoint within each condition relative to time 0
#' @export
#'
#' @examples
mspms_log2fc = function(prepared_for_stats){

  # Extracting just the control data
  control_data = prepared_for_stats %>%
    dplyr::filter(time == 0) %>%
    dplyr::group_by(condition,Peptide) %>%
    dplyr::summarise(control_mean = mean(value,na.rm = TRUE))


  # Extracting just the reference data
  reference_data = prepared_for_stats %>%
    dplyr::group_by(Peptide,time,condition) %>%
    dplyr::summarise(reference_mean = mean(value,na.rm = TRUE))

  # Calculating the log2fc

  log2fc = dplyr::inner_join(control_data,reference_data,by = c("Peptide","condition")) %>%
    dplyr::mutate(comparison = paste0(condition,".T0","_",condition,".T",time),
                  log2fc = log2(reference_mean/control_mean)) %>%
    tibble::as.tibble()


  return(log2fc)


}
