#' mspms_t_tests
#'
#' Performs ttests for each peptide within each group with t0 as reference. P value adjustment performed using FDR correction
#'
#' @param mspms_data = this is the data that has been run through the mspms workflow.
#'
#' @return a tibble with the t test statistics for each peptide within each group with t0 as reference
#' @export
#'
#' @examples
#'
#' t_tests = mspms_t_tests(mspms::mspms_data)
mspms_t_tests = function(mspms_data){

  # Doing T tests
  stat = mspms_data %>%
    dplyr::group_by(.data$cleavage_seq,.data$condition,.data$Peptide) %>%
    rstatix::t_test(value ~ time,ref.group = "0") %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    dplyr::mutate(comparison = paste0(.data$condition,".T",.data$group2,"/",.data$condition,".T0")) %>%
    tibble::as_tibble()

  return(stat)

}
