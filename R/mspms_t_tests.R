#' mspms_t_tests
#'
#' Performs ttests for each peptide within each group with t0 as reference. P value adjustment performed using FDR correction
#'
#' @param prepared_for_stats = this is the data that has been prepared for stats using prepare_for_stats()
#'
#' @return a data frame with the t test statistics for each peptide within each group with t0 as reference
#' @export
#'
#' @examples
#'
#' t_tests = mspms_t_tests(mspms::prepared_for_stats)
mspms_t_tests = function(prepared_for_stats){

  # Doing T tests
  stat = prepared_for_stats %>%
    dplyr::group_by(condition,Peptide) %>%
    rstatix::t_test(value ~ time,ref.group = "0") %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    dplyr::mutate(comparison = paste0(condition,".T0","_",condition,".T",group2)) %>%
    tibble::as_tibble()

  return(stat)

}
