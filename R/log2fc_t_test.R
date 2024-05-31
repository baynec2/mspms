#' log2fc_t_test
#'
#' Gives log2fc and statistics from t tests across conditions. Gives differences relative to time 0 within each condition.
#'
#' @param mspms_data = this is data that has been run through the mspms pipeline.
#'
#' @return a tibble containing log2fc and t test statistics across conditions
#' @export
#'
#' @examples
#'
#' log2fc_and_t_test = log2fc_t_test(mspms::mspms_data)
#'
#'
log2fc_t_test = function(mspms_data){

  # Calculating the log2fc
  log2fc = mspms::mspms_log2fc(mspms_data)

  # Calculating the stats
  stat = mspms::mspms_t_tests(mspms_data)

  #Combining them together
  log2fc_stats = dplyr::inner_join(log2fc,stat,by = c("Peptide","comparison")) %>%
    dplyr::mutate(time = .data$group2,
           condition = .data$condition.x,
           group2 = forcats::fct_inseq(.data$group2)) %>%
    dplyr::select(-.data$condition.x,-.data$condition.y) %>%
    tibble::as_tibble()

  return(log2fc_stats)

}
